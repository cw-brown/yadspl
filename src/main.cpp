#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <typeinfo>
#include <memory>
#include <cassert>

#include "implot.h"

#include "helper_funcs.h"
#include "fir_filter.hpp"
#include "polynomial.hpp"

void key_call(GLFWwindow* window, int key, int, int action, int){
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
}

int main(){
    GLFWwindow* window = glfw_makeNewWindow(1920, 1080, "Yet Another DSP Library", true, true, true);
    ImPlot::CreateContext();
    glfwSetKeyCallback(window, key_call);

    fir_filter filt1, filt2;
    filt1.design_rrc(4, 0.75, 16);
    filt2.design_rc(4, 0.75, 16);

    while(!glfwWindowShouldClose(window)){
        glfwPollEvents();
        glfw_frame();

        ImGui::Begin("Root Raised Cosine Filter");
            if(ImPlot::BeginPlot("Filter")){
                ImPlot::PlotBars("RRC", filt1.getTaps(), filt1.getNumTaps());
                ImPlot::PlotBars("RC", filt2.getTaps(), filt2.getNumTaps());
            ImPlot::EndPlot();
        }
        ImGui::End();

        glfw_render(window);
    }

    glfw_cleanup(window);
    ImPlot::DestroyContext();

    // fir_filter filt;
    // filt.design_rc(4, 0.75, 8);

    // for(size_t i = 0; i < filt.getNumTaps(); ++i)
    //     std::cout<<*(filt.getTaps() + i)<<" ";

    return 0;
}