
/**
 * Some definitions from PortAudio:
 * ================================
 *
 * frame: a frame is a set of samples that occur simultaneously;
 *        for a stereo stream, a frame is two samples.
 *
 * frames per buffer: the number of frames passed to the stream 
 *                    callback function, or the preferred block 
 *                    granularity for a blocking read/write stream.
 */

#include <cmath>
#include <complex>
#include <deque>
#include <iostream> 
#include <list>
#include <memory>
#include <mutex>
#include <sstream>  
#include <string>  
#include <tuple>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include <GLFW/glfw3.h>
#include <fftw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "portaudio.h"




#define PA_SAMPLE_TYPE paFloat32
typedef float SAMPLE_T;
#define SAMPLE_SILENCE  0.0f
#define PRINTF_S_FORMAT "%.8f"

#define NUM_CHANNELS 2 
#define NUM_BUFFERS 1024
#define SAMPLE_RATE  44100
#define FRAMES_PER_BUFFER 2048 
#define SPECTRUM_SIZE 2048

template <typename T>
class SamplePump
{
public:
  
  using BufferT = std::unique_ptr<std::vector<T>>;
  
  SamplePump(
      std::size_t num_frames_per_buffer,
      std::size_t num_channels,
      std::size_t num_max_buffers
  )
  {
    m_num_frames_per_buffer = num_frames_per_buffer;
    m_num_channels = num_channels;

    // Create empty buffers.
    for (std::size_t buff_count=0; buff_count<num_max_buffers; buff_count++)
    {
      m_unused_buffers.push_back(
        std::make_unique<std::vector<T>>(
          m_num_frames_per_buffer * m_num_channels
        )
      );
    }
  }

  BufferT get_empty()
  {
    std::scoped_lock(m_mutex);
    if (m_unused_buffers.empty())
    {
      return nullptr;
    }
    else
    {
      auto r = std::move(m_unused_buffers.back());
      m_unused_buffers.pop_back();
      return r;
    }   
  }

  void return_empty(BufferT buffer)
  {
    std::scoped_lock(m_mutex);
    m_unused_buffers.push_front(std::move(buffer));
  }

  BufferT consume()
  {
    std::scoped_lock(m_mutex);
    if (m_buffers.empty())
    {
      return nullptr;
    }
    else
    {
      auto r = std::move(m_buffers.front());
      m_buffers.pop_front();
      return r;
    }   
  }

  void produce(BufferT buffer)
  {
    std::scoped_lock(m_mutex);
    m_buffers.push_back(std::move(buffer));
  }


  void print_status() const
  {
    printf("%lu unused, %lu used\r",
      m_unused_buffers.size(),
      m_buffers.size()
    );
  }

  std::string get_status() const
  {
    std::stringstream buffer;
    buffer << 
      m_unused_buffers.size() << " unused " <<
      m_buffers.size() << "used";
    return buffer.str(); 
  }

  std::size_t get_num_unused() const
  {
    return m_unused_buffers.size();
  }

private: 
  std::size_t m_num_frames_per_buffer;
  std::size_t m_num_channels;

  std::list<BufferT> m_unused_buffers;
  std::deque<BufferT> m_buffers;
  std::mutex m_mutex;

};


void test_sample_pump()
{
  SamplePump<SAMPLE_T> sp(FRAMES_PER_BUFFER, NUM_CHANNELS, NUM_BUFFERS);
  auto empty = sp.get_empty();
  (*empty)[0] = 5.0;
  (*empty)[1] = 6.0;
  sp.produce(std::move(empty));
  auto full = sp.consume();
  printf("expecting 5.0 and 6.0, got: %1.1f and %1.1f\n\r", (*full)[0], (*full)[1]);
  sp.return_empty(std::move(full));
}


static int record_callback(
  const void* input_buffer, 
  void* output_buffer,
  unsigned long frames_per_buffer,
  const PaStreamCallbackTimeInfo* time_info,
  PaStreamCallbackFlags status_flags,
  void* user_data)
{
    auto input_buffer_ptr = static_cast<const SAMPLE_T*>(input_buffer);
    auto sp = static_cast<SamplePump<SAMPLE_T>*>(user_data);
    auto user_buffer = sp->get_empty();
    if (user_buffer == nullptr)
    {
      // No buffer to process this new data; just skip it.
      return paContinue;
    } 
    auto user_buffer_ptr = user_buffer->data();
    if (input_buffer == nullptr)
    {
      for (unsigned long i=0; i<frames_per_buffer; i++)
      {
        *user_buffer_ptr++ = SAMPLE_SILENCE;
        if (NUM_CHANNELS == 2) 
        {
          *user_buffer_ptr++ = SAMPLE_SILENCE;
        }
      }
    }
    else
    {
      for (unsigned long i=0; i<frames_per_buffer; i++)
      {
        *user_buffer_ptr++ = *input_buffer_ptr++;
        if (NUM_CHANNELS == 2) 
        {
          *user_buffer_ptr++ = *input_buffer_ptr++;
        }
      }
    }

    sp->produce(std::move(user_buffer));
    return paContinue;
}


PaStream* start_portaudio_stream(SamplePump<SAMPLE_T>& sp)
{
  auto err = Pa_Initialize();
  if (err != paNoError)
  {
    throw std::system_error(
      err, 
      std::generic_category(),
      "Could not initialize portaudio"
    );
  }

  PaStream* stream;
  PaStreamParameters inputParameters;
  inputParameters.device = Pa_GetDefaultInputDevice();
  if (inputParameters.device == paNoDevice) 
  {
     throw std::system_error(
      err, 
      std::generic_category(),
      "No default input device"
    );
  }
  inputParameters.channelCount = NUM_CHANNELS;
  inputParameters.sampleFormat = PA_SAMPLE_TYPE;
  inputParameters.suggestedLatency = Pa_GetDeviceInfo(
      inputParameters.device
  )->defaultLowInputLatency;
  inputParameters.hostApiSpecificStreamInfo = NULL;
  
  // Start the pump
  err = Pa_OpenStream(
    &stream,
    &inputParameters,
    NULL, // outputParameters
    SAMPLE_RATE,
    FRAMES_PER_BUFFER,
    paClipOff,      /* we won't output out of range samples so don't bother clipping them */
    record_callback,
    &sp
  );
  
  if (err != paNoError)
  {
    throw std::system_error(
      err, 
      std::generic_category(),
      "Could not open recording stream"
    );
  }

  err = Pa_StartStream(stream);
  if (err != paNoError)
  {
    throw std::system_error(
      err, 
      std::generic_category(),
      "Could not start stream"
    );
  }
  
  return stream;
}


void end_portaudio_stream(PaStream* stream)
{
  PaError err = Pa_CloseStream(stream);
  if (err != paNoError) 
  {
    throw std::system_error(
      err, 
      std::generic_category(),
      "Error while closing stream"
    );
  }
}


static void glfw_error_callback(int error, const char* description)
{
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}


GLFWwindow* start_glfw_and_imgui_context(
  unsigned const int window_size_x, 
  unsigned const int window_size_y
)
{
  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit())
  {
    throw std::system_error(
      -1, 
      std::generic_category(),
      "Error calling glfwInit"
    );
  }

#if defined(IMGUI_IMPL_OPENGL_ES2)
  // GL ES 2.0 + GLSL 100
  const char* glsl_version = "#version 100";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
  // GL 3.2 + GLSL 150
  const char* glsl_version = "#version 150";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
  // GL 3.0 + GLSL 130
  const char* glsl_version = "#version 130";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
  //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

  // Create window with graphics context
  GLFWwindow* window = glfwCreateWindow(
    window_size_x, window_size_y, "Simi", NULL, NULL
  );
  
  if (window == NULL)
  {
    throw std::system_error(
      -1, 
      std::generic_category(),
      "Error calling glfwInit"
    );
  }

  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO(); (void)io;
  //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
  //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  //ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  return window;
}


void end_glfw_and_imgui_context(GLFWwindow* window)
{
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();
}


std::tuple<SAMPLE_T, SAMPLE_T> 
compute_min_max(SamplePump<SAMPLE_T>& sp)
{
  auto buffer = sp.consume();
  if (buffer == nullptr)
  {
    return std::make_tuple(0.0f, 0.0f);
  }
  auto buffer_ptr = buffer->data();
  SAMPLE_T max = *buffer_ptr++;
  SAMPLE_T min = max;
  for (std::size_t i = 1; i<buffer->size(); i++)
  {
    auto v = *buffer_ptr++;
    max = (v > max) ? v : max;
    min = (v < min) ? v : min;
  } 
  sp.return_empty(std::move(buffer));
  return std::make_tuple(min, max);
}


SAMPLE_T global_min = 0.0;
SAMPLE_T global_max = 0.0;


#define V2 ImVec2
void FX(ImDrawList* d, V2 a, V2 b, V2 sz, ImVec4, float t, std::vector<float>& power_spectrum)
{
  float sx = 1.f / (float)power_spectrum.size();
  float sy = 1.f / 9.f;
  int ps_ctr = 0;
  for (float ty = 0.0f; ty < 1.0f; ty += sy)
  {
    for (float tx = 0.0f; tx < 1.0f; tx += sx)
    {
      V2 c((tx + 0.5f * sx), (ty + 0.5f * sy));
      float k = 0.45f;
      d->AddRectFilled(
        V2(a.x + (c.x - k * sx) * sz.x, a.y + (c.y - k * sy) * sz.y),
        V2(a.x + (c.x + k * sx) * sz.x, a.y + (c.y + k * sy) * sz.y),
        IM_COL32(200, 255, 100, power_spectrum[ps_ctr] > ty ? 255 : 50)
      );
    }
    ps_ctr++; 
  }
}


std::vector<SAMPLE_T> average_buffer(FRAMES_PER_BUFFER);
std::vector<std::complex<SAMPLE_T>> fft_buffer(FRAMES_PER_BUFFER);
std::vector<float> power_spectrum(SPECTRUM_SIZE);

void run_imgui_loop(
  GLFWwindow* window, 
  SamplePump<SAMPLE_T>& sp, 
  fftwf_plan& fft_plan
)
{
  // Average left and right into a single channel
  auto buffer = sp.consume();
  
  if (buffer == nullptr)
  {
    // TODO: re-draw window to ensure responsiveness
    return;
  }
  
  auto buffer_ptr = buffer->data();
  auto average_buffer_ptr = average_buffer.data();

  for (std::size_t i = 0; i<FRAMES_PER_BUFFER; i++)
  {
    *average_buffer_ptr++ = (*buffer_ptr++ + *buffer_ptr++) / 2;
  }
  sp.return_empty(std::move(buffer));
  
  // Apply a hanning window to signal
   
  // Compute FFT
  fftwf_execute(fft_plan);

  // compute mangitude: sqrt(real*real + imag*imag)
  // and
  // compute log scale: 20*log10(magnitude)
  // and
  // combine into 8 bins
  
  {
    auto spectrum_ptr = power_spectrum.data();
    auto fft_buffer_ptr = fft_buffer.data();
    auto frames_per_spectrum_bin = (std::size_t)(FRAMES_PER_BUFFER/SPECTRUM_SIZE);
    for (int ps=0; ps<SPECTRUM_SIZE; ps++)
    {
      *spectrum_ptr = 0.0f;
      for (int fq=0; fq<frames_per_spectrum_bin; fq++)
      {
        *spectrum_ptr += 20.0f * std::log10(std::abs(*fft_buffer_ptr++));
      }
      *spectrum_ptr /= (float)frames_per_spectrum_bin; 
      spectrum_ptr++; 
    }   
  }

  
  glfwPollEvents();
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame(); 
 
#if 0 
  ImGui::Begin("Hello, world!");
  ImGui::Text(sp.get_status().c_str());
  ImGui::Text("min %f global min %f", min, global_min);
  ImGui::Text("max %f global max %f", max, global_max);
  ImGui::End();
#elif 0
  ImGuiIO& io = ImGui::GetIO();
  ImGui::Begin("FX", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  ImVec2 size(320.0f, 180.0f);
  ImGui::InvisibleButton("canvas", size);
  ImVec2 p0 = ImGui::GetItemRectMin();
  ImVec2 p1 = ImGui::GetItemRectMax();
  ImDrawList* draw_list = ImGui::GetWindowDrawList();
  draw_list->PushClipRect(p0, p1);

  ImVec4 mouse_data;
  mouse_data.x = (io.MousePos.x - p0.x) / size.x;
  mouse_data.y = (io.MousePos.y - p0.y) / size.y;
  mouse_data.z = io.MouseDownDuration[0];
  mouse_data.w = io.MouseDownDuration[1];

  FX(draw_list, p0, p1, size, mouse_data, (float)ImGui::GetTime(), power_spectrum);
  draw_list->PopClipRect();
  ImGui::End();
#else
  ImGui::Begin("Power Spectrum");
  ImGui::PlotLines(
      "s", 
      power_spectrum.data(), 
      power_spectrum.size(),
      0, // offset
      nullptr, // overlay_text
      -50.0f, // min
      50.0f, // max
      ImVec2(640, 480) // plot size
  ); 
  ImGui::End();
  ImGui::Begin("Input Signal");
  ImGui::PlotLines(
      "i", 
      average_buffer.data(), 
      average_buffer.size(),
      0, // offset
      nullptr, // overlay_text
      -1.0f, // min
      1.0f, // max
      ImVec2(640, 480) // plot size
  ); 
  ImGui::End();
  ImGui::Begin("FFT");
  ImGui::PlotLines(
      "fft", 
      reinterpret_cast<float*>(fft_buffer.data()), 
      fft_buffer.size()*2,
      0, // offset
      nullptr, // overlay_text
      -100.0f, // min
      100.0f, // max
      ImVec2(640, 480) // plot size
  ); 
  ImGui::End();
#endif  

  ImGui::Render();
  int display_w, display_h;
  glfwGetFramebufferSize(window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  glClearColor(
    clear_color.x * clear_color.w, 
    clear_color.y * clear_color.w, 
    clear_color.z * clear_color.w, 
    clear_color.w
  );
  glClear(GL_COLOR_BUFFER_BIT);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  glfwSwapBuffers(window);
}


int main(void)
{
  printf("Simiolus 🐒 is an extinct genus of primates.\n\r");
  test_sample_pump();
  SamplePump<SAMPLE_T> sp(FRAMES_PER_BUFFER, NUM_CHANNELS, NUM_BUFFERS);

  auto window = start_glfw_and_imgui_context(120, 80);

  PaError err = paNoError;
  auto stream = start_portaudio_stream(sp);

  auto fft_plan = fftwf_plan_dft_r2c_1d(
     FRAMES_PER_BUFFER, 
     average_buffer.data(), 
     reinterpret_cast<fftwf_complex*>(fft_buffer.data()),
     FFTW_ESTIMATE
  );

  while ((err = Pa_IsStreamActive(stream)) == 1)
  {
      run_imgui_loop(window, sp, fft_plan);
      // Pa_Sleep(100);
      // sp.print_status();
      // fflush(stdout);
  }

  if (err < 0) 
  {
     throw std::system_error(
      err, 
      std::generic_category(),
      "Error while pumping"
    ); 
  }

  fftwf_destroy_plan(fft_plan); 
  end_portaudio_stream(stream);
  end_glfw_and_imgui_context(window);
  return 0;
}


# if 0


/** @file paex_record.c
    @ingroup examples_src
    @brief Record input into an array; Save array to a file; Playback recorded data.
    @author Phil Burk  http://www.softsynth.com
*/
/*
 * $Id$
 *
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * The text above constitutes the entire PortAudio license; however,
 * the PortAudio community also makes the following non-binding requests:
 *
 * Any person wishing to distribute modifications to the Software is
 * requested to send the modifications to the original developer so that
 * they can be incorporated into the canonical version. It is also
 * requested that these non-binding requests be included along with the
 * license above.
 */

#include <stdio.h>
#include <stdlib.h>
#include "portaudio.h"

/* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
#define SAMPLE_RATE  (44100)
#define FRAMES_PER_BUFFER (512)
#define NUM_SECONDS     (5)
#define NUM_CHANNELS    (2)
/* #define DITHER_FLAG     (paDitherOff) */
#define DITHER_FLAG     (0) /**/
/** Set to 1 if you want to capture the recording to a file. */
#define WRITE_TO_FILE   (0)

/* Select sample format. */
#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif

typedef struct
{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
}
paTestData;

/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int recordCallback( const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData )
{
    paTestData *data = (paTestData*)userData;
    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) outputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;

    if( framesLeft < framesPerBuffer )
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }

    if( inputBuffer == NULL )
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = SAMPLE_SILENCE;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
        }
    }
    else
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}

/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int playCallback( const void *inputBuffer, void *outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo* timeInfo,
                         PaStreamCallbackFlags statusFlags,
                         void *userData )
{
    paTestData *data = (paTestData*)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE*)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) inputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;

    if( framesLeft < framesPerBuffer )
    {
        /* final buffer... */
        for( i=0; i<framesLeft; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
        for( ; i<framesPerBuffer; i++ )
        {
            *wptr++ = 0;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = 0;  /* right */
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else
    {
        for( i=0; i<framesPerBuffer; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}

/*******************************************************************/
int main(void);
int main(void)
{
    PaStreamParameters  inputParameters,
                        outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    int                 i;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;
    SAMPLE              max, val;
    double              average;

    printf("patest_record.c\n"); fflush(stdout);

    data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; /* Record for a few seconds. */
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;
    numBytes = numSamples * sizeof(SAMPLE);
    data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
    if( data.recordedSamples == NULL )
    {
        printf("Could not allocate record array.\n");
        goto done;
    }
    for( i=0; i<numSamples; i++ ) data.recordedSamples[i] = 0;

    err = Pa_Initialize();
    if( err != paNoError ) goto done;

    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto done;
    }
    inputParameters.channelCount = 2;                    /* stereo input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    /* Record some audio. -------------------------------------------- */
    err = Pa_OpenStream(
              &stream,
              &inputParameters,
              NULL,                  /* &outputParameters, */
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              recordCallback,
              &data );
    if( err != paNoError ) goto done;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto done;
    printf("\n=== Now recording!! Please speak into the microphone. ===\n"); fflush(stdout);

    while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
    {
        Pa_Sleep(1000);
        printf("index = %d\n", data.frameIndex ); fflush(stdout);
    }
    if( err < 0 ) goto done;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto done;

    /* Measure maximum peak amplitude. */
    max = 0;
    average = 0.0;
    for( i=0; i<numSamples; i++ )
    {
        val = data.recordedSamples[i];
        if( val < 0 ) val = -val; /* ABS */
        if( val > max )
        {
            max = val;
        }
        average += val;
    }

    average = average / (double)numSamples;

    printf("sample max amplitude = "PRINTF_S_FORMAT"\n", max );
    printf("sample average = %lf\n", average );

    /* Write recorded data to a file. */
#if WRITE_TO_FILE
    {
        FILE  *fid;
        fid = fopen("recorded.raw", "wb");
        if( fid == NULL )
        {
            printf("Could not open file.");
        }
        else
        {
            fwrite( data.recordedSamples, NUM_CHANNELS * sizeof(SAMPLE), totalFrames, fid );
            fclose( fid );
            printf("Wrote data to 'recorded.raw'\n");
        }
    }
#endif

    /* Playback recorded data.  -------------------------------------------- */
    data.frameIndex = 0;

    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto done;
    }
    outputParameters.channelCount = 2;                     /* stereo output */
    outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    printf("\n=== Now playing back. ===\n"); fflush(stdout);
    err = Pa_OpenStream(
              &stream,
              NULL, /* no input */
              &outputParameters,
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              paClipOff,      /* we won't output out of range samples so don't bother clipping them */
              playCallback,
              &data );
    if( err != paNoError ) goto done;

    if( stream )
    {
        err = Pa_StartStream( stream );
        if( err != paNoError ) goto done;

        printf("Waiting for playback to finish.\n"); fflush(stdout);

        while( ( err = Pa_IsStreamActive( stream ) ) == 1 ) Pa_Sleep(100);
        if( err < 0 ) goto done;

        err = Pa_CloseStream( stream );
        if( err != paNoError ) goto done;

        printf("Done.\n"); fflush(stdout);
    }

done:
    Pa_Terminate();
    if( data.recordedSamples )       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    if( err != paNoError )
    {
        fprintf( stderr, "An error occurred while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          /* Always return 0 or 1, but no other return codes. */
    }
    return err;
}

#endif
