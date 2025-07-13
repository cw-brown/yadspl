#include "helper_funcs.h"

#include <exception>
#include <cmath>

GLFWwindow* glfw_makeNewWindow(int width, int height, const char* name, bool fullscreen, bool darkmode, bool vsync){
    if(!glfwInit()) return nullptr;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWmonitor* monitor;
    const GLFWvidmode* mode;
    if(fullscreen){
        monitor = glfwGetPrimaryMonitor();
        mode = glfwGetVideoMode(monitor);
        width = mode->width;
        height = mode->height;
    } else monitor = nullptr;


    GLFWwindow* window = glfwCreateWindow(width, height, name, monitor, nullptr);
    if(!window) return nullptr;

    glfwMakeContextCurrent(window);
    if(vsync) glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    if(darkmode) ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init();

    return window;
}

void glfw_cleanup(GLFWwindow* window){
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}

void glfw_render(GLFWwindow* window){
    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.1f, 0.15f, 0.2f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
}

void glfw_frame(){
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}
