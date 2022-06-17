// ImGui - standalone example application for SDL2 + OpenGL
// If you are new to ImGui, see examples/README.txt and documentation at the top of imgui.cpp.

#include <imgui.h>
#include "imgui_impl_sdl.h"
#include <stdio.h>
#include <cmath>

//#include <thread>
//#include <functional>
//#include <string>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "TinyFileDialogs/tinyfiledialogs.h"

#include"MBdeepzoom/mbimage.h"
#include"MBdeepzoom/mb.h"
#include"IMG2SCR.h"

char const *GetKallesOpenFileName(){
    char const * lTheOpenFileName;
    char const * lFilterPatterns[2] = { "*.kfr" };
    lTheOpenFileName = tinyfd_openFileDialog(
            "Open kalles fraktaler file",
            "",
            1,
            lFilterPatterns,
            NULL,
            0);
    return lTheOpenFileName;
}

char const *GetKallesSaveFileName(){
    char const * lTheSaveFileName;
    char const * lFilterPatterns[2] = { "*.kfr" };
    lTheSaveFileName = tinyfd_saveFileDialog(
            "Save kalles fraktaler file",
            "",
            1,
            lFilterPatterns,
            NULL);
    return lTheSaveFileName;
}

char const *GetPPMSaveFileName(){
    char const * lTheSaveFileName;
    char const * lFilterPatterns[2] = { "*.ppm" };
    lTheSaveFileName = tinyfd_saveFileDialog(
            "Save to ppm",
            "",
            1,
            lFilterPatterns,
            NULL);
    return lTheSaveFileName;
}

int SDLCALL doSomething(void *data){
    //Uint32 t0 = SDL_GetTicks();//looks like SDL_GetTicks() is not thread safe!
    MB *pMB = (MB *) data;
    pMB->dowork();
    //pMB->m_time = SDL_GetTicks() - t0;
   return 0;
}

void textureFromImage(MBImage &img){
    void *ptr = img.getPtr2Img();
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width,
    img.height, 0, GL_RGB, GL_UNSIGNED_BYTE,
    ptr);
}

void drawScrSimple(ImVec4 &clear_color, GLuint texName, MBImage &img);
void drawScr(ImVec4 &clear_color, GLuint texName, ImgScr &CConverter);
void drawFrames(float mx, float my, float zoomBy, ImgScr &CConverter);

int main(int, char**)

{
    //zero location
    const char *sre = "0.";
    const char *sim = "0.";
    long double sradius = 2.E0L;
    int precision = 70;
    int maxIter   = 500;
    int seriesApproxOrder = 32;

    int skip = -1;
    float maxSAerr = 1.;
    float maxGCerr = 0.001;
    bool solveGlitches = true;
    int GD_method = GD_GERRIT_B;
    const char* GD_text[6] = {"None","Pdb","Kn","KnPdb","Gerrit","GerritB"};
    int maxGlitchPasses = 128;


    int imgWidth = 640, imgHeight = 360;

    MBImage img(imgWidth,imgHeight);

    MB MBinst;
    MBinst.setPrecision(precision);
    MBinst.setLocation(sre,sim);
    MBinst.setSRadius(sradius);
    MBinst.setMaxIter(maxIter);
    MBinst.setSeriesApproxOrder(seriesApproxOrder);
    MBinst.setSkipIter(skip);
    MBinst.m_allowedErrorSA = maxSAerr;
    MBinst.m_allowedErrorGC = maxGCerr;
    MBinst.m_solveGlitches = solveGlitches;
    MBinst.setImage(img);


    // Setup SDL
    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0)
    {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }
    // Setup window
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    SDL_DisplayMode current;
    SDL_GetCurrentDisplayMode(0, &current);
    SDL_Window *window = SDL_CreateWindow("Mandelbrot set deep zoom experiment", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
    //SDL_Surface* screenSurface = NULL;
    SDL_GLContext glcontext = SDL_GL_CreateContext(window);

    // Setup ImGui binding
    ImGui_ImplSdl_Init(window);

    // Load Fonts
    // (there is a default font, this is only if you want to change it. see extra_fonts/README.txt for more details)
    //ImGuiIO& io = ImGui::GetIO();
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("../../extra_fonts/Cousine-Regular.ttf", 15.0f);
    //io.Fonts->AddFontFromFileTTF("../../extra_fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../extra_fonts/ProggyClean.ttf", 13.0f);
    //io.Fonts->AddFontFromFileTTF("../../extra_fonts/ProggyTiny.ttf", 10.0f);
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());

    bool show_test_window = false;
    bool show_another_window = false;
    ImVec4 clear_color = ImColor(114, 144, 154);

    //OpenGL init
    GLuint texName;//We have only one texture
    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);//GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);//GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    textureFromImage(img);

    // Main loop
    bool done = false;
    Uint32 t0=0,t1=1;//used to limit FPS
    Uint32 rt0=0, rt1=1, rt2=0; //used for rendering time statistics
    bool reloadImg = false;
    bool bShowRoots = false;
    bool bDoRender = true;
    /*static*/ int rootnum=0;
    float iZoomBy = 2.f;
    bool snap2root = false;
    //Uint32 prevImgNbr = 0;
    ImgScr CConverter(ImGui::GetIO().DisplaySize.x,ImGui::GetIO().DisplaySize.y);
    CConverter.SetTexDims(img.width,img.height);
    SDL_ShowCursor(0);
    while (!done)
    {
        t0 = SDL_GetTicks();
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            ImGui_ImplSdl_ProcessEvent(&event);
            if (event.type == SDL_QUIT)
                done = true;
        }

        ImGui_ImplSdl_NewFrame(window);
        ImGui::GetIO().MouseDrawCursor=true;

        CConverter.SetScreenDims(ImGui::GetIO().DisplaySize.x,ImGui::GetIO().DisplaySize.y);
        //CConverter.SetTcenPos(0.5f*ImGui::GetIO().DisplaySize.x,0.5f*ImGui::GetIO().DisplaySize.y);
        CConverter.Tcenter();

        // 1. Show a simple window
        // Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appears in a window automatically called "Debug"
        {
            ImGui::SetNextWindowSize(ImVec2(200,100), ImGuiSetCond_FirstUseEver);
            ImGui::Begin("Mandelbrot set deep zoom experiment", &show_another_window);
            //ImGui::Text("Mandelbrot set deep zoom experiment!");
            if(MBinst.m_status == MB::S_idle){
                if (ImGui::Button("Zoom in")) {
                    sradius *= 1./iZoomBy;//.5;
                    bDoRender = true;
                }
                ImGui::SameLine();
                if (ImGui::Button("Zoom out")) {
                    sradius *= iZoomBy;//2.;
                    bDoRender = true;
                }
                if (ImGui::Button("Render")){
                    bDoRender = true;
                }
            } else {
                if (ImGui::Button("Abort")) {
                    MBinst.m_request = MB::R_abort;
                }
            }

            if(MBinst.m_status == MB::S_idle){
                if(ImGui::CollapsingHeader("File")){
                    if(ImGui::Button("Home")){
                        sradius = 2.E0L;
                        precision = 70;
                        maxIter   = 500;
                        seriesApproxOrder = 32;
                        imgWidth = 640; imgHeight = 360;
                        MBinst.setPrecision(precision);
                        MBinst.setLocation(sre,sim);
                        MBinst.setSRadius(sradius);
                        MBinst.setMaxIter(maxIter);
                        MBinst.setSeriesApproxOrder(seriesApproxOrder);
                        MBinst.setSkipIter(-1);
                        bDoRender = true;
                    }
                    if(ImGui::Button("Read kalles fractaler")){
                        const char *fName = GetKallesOpenFileName();
                        if(fName)
                            if(MBinst.load(fName)){
                                sradius = MBinst.getSRadius();
                                maxIter = MBinst.getMaxIter();
                                skip = -1;
                                seriesApproxOrder = 64;
                                SDL_SetWindowTitle(window,fName);
                            }
                    }
                    if(ImGui::Button("Save kalles fractaler")){
                        const char *fName = GetKallesSaveFileName();
                        if(fName)
                            if(MBinst.save(fName)){
                                SDL_SetWindowTitle(window,fName);
                            }
                    }
                }
            }
            if(ImGui::CollapsingHeader("Paramaters")){
                ImGui::InputInt("Skip", &skip);
                ImGui::InputInt("Approx. Order", &seriesApproxOrder); seriesApproxOrder = std::min(std::max(seriesApproxOrder,4),256);
                ImGui::InputInt("Max. Iter.", &maxIter); maxIter = std::max(100, maxIter);
            }

            if(ImGui::CollapsingHeader("Navigation")){
                ImGui::InputFloat("Zoom by", &iZoomBy); iZoomBy = std::min(16.f,std::max(1.f,iZoomBy));
                {
                    float mantissa;
                    int   expo;
                    mantissa = frexpl(sradius,&expo);
                    ImGui::SliderFloat("mant.", &mantissa, 0.5, 1.0);
                    ImGui::InputInt("2^expo.", &expo);
                    sradius = ldexpl(mantissa,mantissa == 1.f ? expo-1 : expo);
                }
                ImGui::Checkbox("Show roots",&bShowRoots);
                ImGui::SameLine();
                ImGui::Checkbox("Snap to roots",&snap2root);
            }

            if(ImGui::CollapsingHeader("Series accuracy")){
                ImGui::InputFloat("Max SA err", &maxSAerr); maxSAerr = std::min(1000.f, std::max(0.0001f,maxSAerr));
            }
            if(ImGui::CollapsingHeader("Glitch handling")){
                if(ImGui::Combo("Glitch detection method", &GD_method, GD_text, 6)){
                    if(GD_method == GD_KNPDB)
                        maxGCerr = 0.000001;
                    if(GD_method == GD_KN || GD_method == GD_GERRIT)
                        maxGCerr = 0.0001;
                    if(GD_method == GD_GERRIT_B)
                        maxGCerr = 0.001;
                }
                ImGui::Checkbox("Solve Glitches",&solveGlitches);
                ImGui::InputInt("Max GC passes", &maxGlitchPasses); maxGlitchPasses = std::min(std::max(maxGlitchPasses,0),256);
                ImGui::InputFloat("Max Glitch err", &maxGCerr); maxGCerr = std::min(1000.f, std::max(0.00000001f,maxGCerr));
            }
            if(ImGui::CollapsingHeader("Display"))//CollapsingHeader("Display"))
            {
                ImGui::ColorEdit3("clear color", (float*)&clear_color);
                ImGui::InputInt("Image width", &imgWidth); imgWidth = std::min(std::max(imgWidth,64),16384);
                ImGui::InputInt("Image height", &imgHeight); imgHeight = std::min(std::max(imgHeight,64),16384);
                if(ImGui::Button("Save to .PPM")){
                    const char *fName = GetPPMSaveFileName();
                    if(fName) img.save(fName);
                }
            }

            if(ImGui::CollapsingHeader("Information"))
            {
              int expo = int(std::log10(sradius));
              ImGui::Text("Zoom factor : %1.3gE%d", double(sradius * std::pow(10.L, -expo)), expo);
              ImGui::Text("Series iters : %d", MBinst.m_seriesNbr);
              ImGui::Text("Reference iters : %d", MBinst.m_referenceNbr);
              ImGui::Text("Rendered lines : %d", MBinst.m_imageNbr);
              ImGui::Text("Glitches : %d", MBinst.m_glitchNbr);
              ImGui::Text("GC passes : %d", MBinst.m_glitchPasses);
              ImGui::Text("Min iterations : %d", MBinst.m_minIterNbr);
              ImGui::Text("Max iterations : %d", MBinst.m_maxIterNbr);
              int rootsFound = MBinst.getNbrRoots();
              ImGui::Text("Roots found : %d", rootsFound);
#if 0
              ImGui::Text("Time : %.3f s", MBinst.m_time * 0.001f);
#else
              {
                  unsigned int tics,sec, min, hrs;//, days;
                  tics = MBinst.m_time;
                  tics /= 10;
                  sec   = tics / 100; tics = tics % 100;
                  min   = sec  / 60;  sec  = sec  %  60;
                  hrs   = min  / 60;  min  = min  %  60;
                  //days  = hrs  / 24;  hrs  = hrs  %  24;
                  ImGui::Text("Time : %3d:%02d:%02d.%02d ", hrs, min, sec, tics);
              }
#endif
              //ImGui::Text("io.DeltaTime : %g", io.DeltaTime);
              //ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
              //ImGui::Text("Application average %d ms", t1-t0);
              if(rootsFound > 0)
              if(ImGui::TreeNode("Roots :"))//CollapsingHeader("Info."))
              {
                  ImGui::InputInt("Root index", &rootnum); rootnum = std::min(std::max(rootnum,0),rootsFound-1);
                  long double r_x,r_y;
                  int period, quality;
                  MBinst.getRoot(rootnum,r_x,r_y, period, quality);
                  ImGui::Text("Period  : %d", period);

                  expo = r_x==0.L ? 0 : int(std::log10(std::abs(r_x)));
                  ImGui::Text("Rel X   : %1.6gE%d", double(r_x * std::pow(10.L, -expo)), expo);
                  expo = r_y==0.L ? 0 : int(std::log10(std::abs(r_y)));
                  ImGui::Text("Rel Y   : %1.6gE%d", double(r_y * std::pow(10.L, -expo)), expo);

                  ImGui::Text("Quality : %d", quality);
                  ImGui::TreePop();
              }
            }

            if (show_test_window)
                    {
                        ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
                        ImGui::ShowTestWindow(&show_test_window);
                    }

            ImGui::End();
        }

//Set the cursor types
        if(ImGui::GetIO().WantCaptureMouse)
            ImGui::SetMouseCursor(ImGuiMouseCursor_Arrow);
        else
            ImGui::SetMouseCursor(ImGuiMouseCursor_None);

//Clickn'zoom
        CConverter.Tcenter();
        if(MBinst.m_status == MB::S_idle && !ImGui::GetIO().WantCaptureMouse){
            ImVec2 sp = ImGui::GetMousePos();
            sp.y = ImGui::GetIO().DisplaySize.y - sp.y;
            float scl = 1.f;
            if (ImGui::IsMouseReleased(0))
                scl = 1./iZoomBy;
            else if (ImGui::IsMouseReleased(1))
                scl = iZoomBy;
            if(ImGui::IsMouseReleased(0) || ImGui::IsMouseReleased(1)){//scl != 1.f){
                float rx = sp.x, ry = sp.y;
                CConverter.S2T(rx,ry,rx,ry);
                long double lx = rx, ly = ry;
                //use root near cursor if distance to cursor less than 6 pixels
                if(snap2root && MBinst.getNbrRoots()>0){
                    int nri = -1;
                    for(unsigned int i=0; i < MBinst.getNbrRoots(); i++){
                        float rx,ry;
                        MBinst.getRootF(i,rx,ry);
                        CConverter.T2S(rx,ry,rx,ry);
                        float d = std::max(std::abs(rx-sp.x),std::abs(ry-sp.y));
                        if(d < 6){
                            nri  = i;
                            break;
                        }
                    }
                    if (nri >= 0){
                        int dummy, dummy1;
                        MBinst.getRoot(nri,lx,ly,dummy,dummy1);
                        MBinst.win2tex(lx,ly,lx,ly);
                    }
                }
                //---------------------
                sradius *= scl;
                MBinst.setLocRel(lx, ly, scl);
                bDoRender = true;
            }
        }
        //Set zoom factor using mousewheel
        if( !ImGui::GetIO().WantCaptureMouse){
            int m = int(ImGui::GetIO().MouseWheel);
            if(m > 0)
                iZoomBy = std::min(iZoomBy * 1.090507733, 16.);
            else if(m<0)
                iZoomBy = std::max(iZoomBy / 1.090507733, 1.);
        }

//Compute picture
        if (MBinst.m_status == MB::S_idle && bDoRender){
            bDoRender = false;
            //
            img.resize(imgWidth,imgHeight);
            MBinst.setSkipIter(skip);
            MBinst.setMaxIter(maxIter);
            MBinst.setSeriesApproxOrder(seriesApproxOrder);
            MBinst.setSRadius(sradius);

            MBinst.m_allowedErrorSA = maxSAerr;
            MBinst.m_allowedErrorGC = maxGCerr;
            MBinst.m_solveGlitches = solveGlitches;
            MBinst.m_GDmethod = GD_method;
            MBinst.m_maxGlitchPasses = maxGlitchPasses;
            MBinst.m_request = MB::R_doWork;

            rt0 = SDL_GetTicks();
            #if 1
            //Use SDL threads for the MB set rendering: Simple and works perfectly.
            SDL_Thread *thread;
            thread = SDL_CreateThread(doSomething, "MB", &MBinst);
            SDL_DetachThread(thread);
            #else
            //Works well... but obviously it hangs the GUI.
            doSomething(&MBinst);
            #endif
            reloadImg = true;
            //prevImgNbr = 0;
        }

        if(MBinst.m_status == MB::S_busy){
            rt1 = SDL_GetTicks();
            MBinst.m_time = rt1 - rt0;
        }
// Rendering screen
        //Do we have a new render?
        if(reloadImg && (MBinst.m_status == MB::S_finished || MBinst.m_status == MB::S_busy)){
            rt1 = SDL_GetTicks();
            if(MBinst.m_status == MB::S_finished){
                reloadImg = false;

                MBinst.m_time = rt1 - rt0;

                MBinst.m_status = MB::S_idle;
            }
            if(MBinst.m_status == MB::S_busy && /*(MBinst.m_imageNbr - prevImgNbr < 20 && MBinst.m_imageNbr < img.height) */(int(rt1) - int(rt2)) < 500){

            } else {
                glBindTexture(GL_TEXTURE_2D, texName);
                //if(bShowRoots) MBinst.drawRoots();
                textureFromImage(img);
                CConverter.SetTexDims(img.width,img.height);
                rt2 = rt1;
                //prevImgNbr = MBinst.m_imageNbr;
            }
        }
        //CConverter.SetTcenPos(0.5f*ImGui::GetIO().DisplaySize.x,0.5f*ImGui::GetIO().DisplaySize.y);
        ////CConverter.Tcenter();
        //drawScrSimple(clear_color, texName, img);
        drawScr(clear_color, texName, CConverter);
#if 1
        //Draw roots
        if(bShowRoots){
            float xp,yp;
            float dx = 2./ImGui::GetIO().DisplaySize.x * 4.;
            float dy = 2./ImGui::GetIO().DisplaySize.y * 4.;
            glColor3f(0.,0.5,1.);
            if(MBinst.getNbrRoots()>0){
                int i = 0;
                while(MBinst.getRootF(i,xp,yp)){
                    if(i != rootnum){
                        CConverter.T2S1(xp,yp,xp,yp);
                        glBegin(GL_LINE_LOOP);
                            glVertex2f(xp-dx, yp-dy);
                            glVertex2f(xp+dx, yp+dy);
                            glVertex2f(xp-dx, yp+dy);
                            glVertex2f(xp+dx, yp-dy);
                        glEnd();
                    }
                    i++;
                }
                glColor3f(1.,0.25,0.);
                MBinst.getRootF(rootnum,xp,yp);
                CConverter.T2S1(xp,yp,xp,yp);
                glBegin(GL_LINE_LOOP);
                    glVertex2f(xp-dx, yp-dy);
                    glVertex2f(xp+dx, yp+dy);
                    glVertex2f(xp-dx, yp+dy);
                    glVertex2f(xp+dx, yp-dy);
                glEnd();
            }
        }

        //Draw zooming frame
        if(!ImGui::GetIO().WantCaptureMouse)
            drawFrames(ImGui::GetIO().MousePos.x, ImGui::GetIO().MousePos.y, iZoomBy, CConverter);

        glColor3f(1.,1.,1.);
#endif
        //Draw GUI
        ImGui::Render();
        SDL_GL_SwapWindow(window);

        //Limit framerate to roughly 50 fps...
        t1 = SDL_GetTicks();
        int delay = 20 - (int(t1) - int(t0));
        if(delay>0) SDL_Delay(Uint32(delay));
    }

    // Cleanup
    // TODO: In case rendering took too long and the user closes the window,
    //       we should tell the working thread to stop and exit.

    ImGui_ImplSdl_Shutdown();
    SDL_GL_DeleteContext(glcontext);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

void drawScrSimple(ImVec4 &clear_color, GLuint texName, MBImage &img){
    glEnable(GL_TEXTURE_2D);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, &clear_color.x);
    glViewport(0, 0, (int)ImGui::GetIO().DisplaySize.x, (int)ImGui::GetIO().DisplaySize.y);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindTexture(GL_TEXTURE_2D, texName);
    glBegin(GL_QUADS);
        float sw = ImGui::GetIO().DisplaySize.x, sh = ImGui::GetIO().DisplaySize.y;
        float cx = floor(sw*.5), cy = floor(sh*.5);//For now
        float Tx0 = 0.5 - cx / img.width, Tx1 = (sw - cx)/img.width + .5;
        float Ty1 = 0.5 - cy / img.height, Ty0 = (sh - cy)/img.height + .5;
        glTexCoord2f(Tx0, Ty0); glVertex2f(-1., -1.);
        glTexCoord2f(Tx1, Ty0); glVertex2f(1., -1.);
        glTexCoord2f(Tx1, Ty1); glVertex2f(1., 1.);
        glTexCoord2f(Tx0, Ty1); glVertex2f(-1., 1.);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

void drawScr(ImVec4 &clear_color, GLuint texName, ImgScr &CConverter){
    float sw = ImGui::GetIO().DisplaySize.x;
    float sh = ImGui::GetIO().DisplaySize.y;
    //float tw(img.width), th(img.height);
    float tw(CConverter.m_tWidth), th(CConverter.m_tHeight);
    glEnable(GL_TEXTURE_2D);
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, &clear_color.x);
    glViewport(0, 0, (int)sw, (int)sh);
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindTexture(GL_TEXTURE_2D, texName);
    glBegin(GL_QUADS);
        float Tx , Ty;
        CConverter.S2T(0,0,Tx,Ty); Tx/=tw; Ty/=th;
        glTexCoord2f(Tx, Ty); glVertex2f(-1., -1.);
        CConverter.S2T(sw,0,Tx,Ty); Tx/=tw; Ty/=th;
        glTexCoord2f(Tx, Ty); glVertex2f(1., -1.);
        CConverter.S2T(sw,sh,Tx,Ty); Tx/=tw; Ty/=th;
        glTexCoord2f(Tx, Ty); glVertex2f(1., 1.);
        CConverter.S2T(0,sh,Tx,Ty); Tx/=tw; Ty/=th;
        glTexCoord2f(Tx, Ty); glVertex2f(-1., 1.);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

void drawFrames(float mx, float my, float zoomBy, ImgScr &CConverter){
    my = CConverter.m_sHeight - my;
    CConverter.S2T(mx,my,mx,my);
    //Draw frame. Badly written. ToDo: a frames hyerarchy that is understandable!!!
    float tw = CConverter.m_tWidth / zoomBy;
    float th = CConverter.m_tHeight / zoomBy;
    float x,y,w,h,pw,ph;
    CConverter.S2T(0.f,0.f,x,y); CConverter.S2T(tw,th,w,h); w-=x; h-=y; w*=0.5; h*=0.5;
    CConverter.S2T(0.f,0.f,x,y); CConverter.S2T(1.f,1.f,pw,ph); pw-=x; ph-=y;

    glColor3f(1.f,1.f,1.f);
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx-w-pw,my-h-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w+pw,my-h-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w+pw,my+h+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-w-pw,my+h+ph,x,y); glVertex2f(x,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx-w+pw,my-h+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w-pw,my-h+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w-pw,my+h-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-w+pw,my+h-ph,x,y); glVertex2f(x,y);
    glEnd();

    glColor3f(0.f,0.f,0.f);
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx-w,my-h,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w,my-h,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w,my+h,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-w,my+h,x,y); glVertex2f(x,y);
    glEnd();

    CConverter.S2T(0.f,0.f,x,y); CConverter.S2T(8.f,8.f,w,h); w-=x; h-=y;
    glColor3f(0.f,0.f,0.f);
    glBegin(GL_LINES);
        CConverter.T2S1(mx+2.f*pw,my,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w,my,x,y); glVertex2f(x,y);

        CConverter.T2S1(mx-2.f*pw,my,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-w,my,x,y); glVertex2f(x,y);

        CConverter.T2S1(mx,my+2.f*pw,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx,my+h,x,y); glVertex2f(x,y);

        CConverter.T2S1(mx,my-2.f*pw,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx,my-h,x,y); glVertex2f(x,y);
    glEnd();

    glColor3f(1.f,1.f,1.f);
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx+pw,my-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w+pw,my-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+w+pw,my+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+pw,my+ph,x,y); glVertex2f(x,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx+pw,my+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+pw,my+h+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-pw,my+h+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-pw,my+ph,x,y); glVertex2f(x,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx-pw,my-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-w-pw,my-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-w-pw,my+ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-pw,my+ph,x,y); glVertex2f(x,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
        CConverter.T2S1(mx+pw,my-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx+pw,my-h-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-pw,my-h-ph,x,y); glVertex2f(x,y);
        CConverter.T2S1(mx-pw,my-ph,x,y); glVertex2f(x,y);
    glEnd();

    glColor3f(1.f,1.f,1.f);
}
