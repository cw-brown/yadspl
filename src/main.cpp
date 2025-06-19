#include <iostream>
#include <complex>
#include <array>
#include <memory>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <utility>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
#include "iir_filter.hpp"
#include "polynomial.hpp"
#include "noise.hpp"
#include "yadpsl_math.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

#define DO_WINDOW false

int main(){
    using namespace std::literals::complex_literals;
    complex_polynomial<double> b = {1.0 + 1i, 2.0 - 1i};
    complex_polynomial<double> c = {3.0 - 3i};

    complex_polynomial<double> a(4);
    a.fill(2.0 + 2.0i);
    std::cout<<a.front();

    std::cout<<"Order: "<<a.order()<<"\n";
    for(auto&& v : a){
        std::cout<<v<<" ";
    }


    if(DO_WINDOW){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    while(!glfwWindowShouldClose(window)){
        glfwPollEvents();
        glfw_frame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, ImGui::GetIO().DisplaySize.y));
        ImGuiWindowFlags topbarflags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse;

        ImGui::Begin("Plottings", nullptr, topbarflags);
        ImGui::BeginTabBar("Main Tabs");
        if(ImGui::BeginTabItem("Test Feature")){
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();
    }
    return 0;
}

