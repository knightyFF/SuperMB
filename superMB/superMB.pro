TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

IM_GUI_INCLUDE_PATH = $$PWD"/imgui"
#message($$IM_GUI_INCLUDE_PATH)

INCLUDEPATH += $$IM_GUI_INCLUDE_PATH

SOURCES += \
    imgui_impl_sdl.cpp \
    main.cpp \
    imgui/imgui.cpp \
    imgui/imgui_demo.cpp \
    imgui/imgui_draw.cpp \
    MBdeepzoom/MBImage.cpp \
    MBdeepzoom/mb.cpp \
    TinyFileDialogs/tinyfiledialogs.c

HEADERS += \
    imgui/imgui.h \
    imgui_impl_sdl.h \
    imgui/imconfig.h \
    imgui/imgui_internal.h \
    imgui/stb_rect_pack.h \
    imgui/stb_textedit.h \
    imgui/stb_truetype.h \
    MBdeepzoom/complex.h \
    MBdeepzoom/mbimage.h \
    MBdeepzoom/mb.h \
    TinyFileDialogs/tinyfiledialogs.h \
    IMG2SCR.h \
    mpfr/mpreal.h

##Minimal library set to link with in order to use SDL2
LIBS += -lmingw32
LIBS += -lSDL2main -lSDL2.dll
LIBS += -lopengl32
##Other libraries
LIBS += -lmpfr -lgmp -fopenmp ##-lgomp ## So... "-fopenmp" needs to be in the linker options too! or add "-lgomp"
win32 {
    ##For tinyFileDialog on windows
    LIBS += -lComdlg32 -lOle32
}

QMAKE_CXXFLAGS += -std=c++11 -Ofast -mfpmath=387 -march=native -fopenmp -Wall -Wextra -pedantic
##QMAKE_CXXFLAGS_RELEASE -= -O2
##QMAKE_CXXFLAGS -= -O2
QMAKE_CXXFLAGS_RELEASE -= -fexceptions
QMAKE_CXXFLAGS_DEBUG -= -Ofast
QMAKE_CXXFLAGS_DEBUG += -O0
##QMAKE_CXXFLAGS_RELEASE += -Ofast -fno-signed-zeros -fno-trapping-math
