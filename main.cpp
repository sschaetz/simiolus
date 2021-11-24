
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

#include <algorithm>
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
#define SPECTRUM_SIZE 512


std::list<std::tuple<int, int>> equalizer =
{
  {5, 10},
  {10, 15},
  {15, 20},
  {20, 25},
  {25, 30},
  {30, 35},
  {35, 40},
  {40, 45},
  {45, 50},
  {50, 55},
  {55, 60},
  {60, 65},
  {65, 70},
  {70, 75},
  {75, 80},
  {80, 85},
  {85, 90},
  {90, 95},
  {95, 100},
  {100, 105},
  {105, 110},
  {110, 115},
  {115, 120},
  {120, 125},
  {125, 130},
  {130, 135},
  {135, 140},
  {140, 145},
  {145, 150},
  {150, 155},
  {155, 160},
  {160, 165},
  {165, 170},
  {170, 175},
  {175, 180},
  {180, 185},
  {185, 190},
  {190, 195},
  {195, 200}
};


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
    auto l = std::scoped_lock(m_mutex);
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
    auto l = std::scoped_lock(m_mutex);
    m_unused_buffers.push_front(std::move(buffer));
  }

  BufferT consume()
  {
    auto l = std::scoped_lock(m_mutex);
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
    auto l = std::scoped_lock(m_mutex);
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
std::vector<SAMPLE_T> hann_filter(FRAMES_PER_BUFFER);
std::vector<std::complex<SAMPLE_T>> fft_buffer(FRAMES_PER_BUFFER);
std::vector<float> power_spectrum(SPECTRUM_SIZE);
std::vector<float> equalizer_buffer(equalizer.size());


void compute_hann_filter(std::vector<SAMPLE_T>& hann_filter)
{
  auto N = hann_filter.size()-1;
  for (std::size_t i=0; i<hann_filter.size(); i++)
  {
    hann_filter[i] = .5 * (1 - std::cos(2*M_PI*i/N));
  }
}


void compute_equalizer()
{
  int start, stop;
  auto equlizer_buffer_ptr = equalizer_buffer.data();
  auto power_spectrum_ptr = power_spectrum.data();
  for(auto range : equalizer)
  {
    std::tie(start, stop) = range;
    *equlizer_buffer_ptr = 0.0f;
    for (int i=start; i<stop; i++)
    {
      *equlizer_buffer_ptr += -1.0f * std::clamp(
          *power_spectrum_ptr++,
          -100.0f,
          0.0f
        );
    }
    *equlizer_buffer_ptr /= (float)(stop-start);
    // printf("eq %d %d %f\n", start, stop, *equlizer_buffer_ptr);
    equlizer_buffer_ptr++;
  }
}


#define V2 ImVec2

void run_imgui_loop(
  GLFWwindow* window,
  SamplePump<SAMPLE_T>& sp,
  fftwf_plan& fft_plan
)
{
  // Average left and right into a single channel
  // and
  // Apply a hanning window to signal
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
    *average_buffer_ptr++ = ((*buffer_ptr++ + *buffer_ptr++) / 2) *
      hann_filter[i];
  }
  sp.return_empty(std::move(buffer));

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
    auto frames_per_spectrum_bin = (std::size_t)(FRAMES_PER_BUFFER/SPECTRUM_SIZE/2);
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

  compute_equalizer();

  glfwPollEvents();
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  ImGui::Begin("Power Spectrum");
  ImGui::PlotLines(
      "s",
      power_spectrum.data(),
      power_spectrum.size(),
      0, // offset
      nullptr, // overlay_text
      -50.0f, // min
      50.0f, // max
      ImVec2(230, 180) // plot size
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
      ImVec2(550, 200) // plot size
  );
  ImGui::End();
  ImGui::Begin("FFT");
  ImGui::PlotLines(
      "fft",
      reinterpret_cast<float*>(fft_buffer.data()),
      fft_buffer.size(),
      0, // offset
      nullptr, // overlay_text
      -100.0f, // min
      100.0f, // max
      ImVec2(550, 200) // plot size
  );
  ImGui::End();

  ImGuiIO& io = ImGui::GetIO();
  ImGui::Begin("Equalizer", NULL, ImGuiWindowFlags_AlwaysAutoResize);
  ImVec2 size(320.0f, 180.0f);
  ImGui::InvisibleButton("canvas", size);
  ImVec2 a = ImGui::GetItemRectMin();
  ImVec2 b = ImGui::GetItemRectMax();
  ImDrawList* draw_list = ImGui::GetWindowDrawList();
  draw_list->PushClipRect(a, b);

  {
    float sx = 1.f / (float)equalizer_buffer.size();
    float sy = 1.f / 16.0f;
    int ps_ctr = 0;
    auto max = *std::max_element(
      equalizer_buffer.begin(),
      equalizer_buffer.end()
    );
    for (float tx = 0.0f; tx < 1.0f; tx += sx)
    {
      // This goes from left to right.
      auto eq_val = 1.0f - (equalizer_buffer[ps_ctr] / max);
      for (float ty = 0.0f; ty < 1.0f; ty += sy)
      // This goes from top to bottom.
      {
        V2 c((tx + 0.4f * sx), (ty + 0.4f * sy));
        float k = 0.3f;
        draw_list->AddRectFilled(
          V2(a.x + (c.x - k * sx) * size.x, a.y + (c.y - k * sy) * size.y),
          V2(a.x + (c.x + k * sx) * size.x, a.y + (c.y + k * sy) * size.y),
          IM_COL32(200, 255, 100, eq_val >= 1.0f - ty ? 255 : 50)
        );
      }
      ps_ctr++;
    }
  }

  draw_list->PopClipRect();
  ImGui::End();

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

  compute_hann_filter(hann_filter);
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

