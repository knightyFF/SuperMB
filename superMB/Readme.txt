SuperMB is a simple somewhat deep mandelbrot set zooming program.
It uses ImGUI and TinyFileDialogs for GUI and mpreal with mpfr and GMP for high precision computation.
By Abdelaziz Naït (knighty) Merzouk. Based on original code by Claude (claude) Heiland-Allen. (2016 - 2018)
No waranty of any kind. Use at your own risk.

3rd party code:
* mpreal: in mpfr folder. Links: See the .h file in that folder.
* imGUI : in imgui folder. main.cpp was developed from one of the usage examples of imGUI. Links: see the readme.md in that folder.
* Tiny File Dialog: in TinyFileDialogs. Links: See the readme in that folder.

Capabilities: 
This is highly experimental. It was first intended to explore ways to improve superfractalthing method: Series approximation and glitch detection.
* Deep zooming limited by long double floats range. The deepest for now was around 1e-2400. In principle up to 1e-4000 (maybe 1e-8000 ?).
* Click to zoom: left mouse button to zoom in, right mouse button to zoom out and mouse wheel to change the zoomin factor per click.
* Loads and saves to kalles fraktaler location file type.
* Saves picture in .ppm format. You need an external application to convert to other image formats.
* Shows some of the periodic points (root). Those are used internally but can also help looking at where to zoom in. Option to snap to a root near cursor while zooming.
* Renders only in distance estimation mode for now.

Compiling:
The author don't know a lot about compiling. Use msys2 + minGW + GCC and Qt creator as IDE:
- install msys2 + minGW + GCC if you don't have them. The guide: https://github.com/orlp/dev-on-windows/wiki/Installing-GCC--&-MSYS2
- install msys2 + minGW + GCC in Qt creator as a new kit. Menu/tools/Options then in the dialog choose "compilation & execution" in the left widget... then find your way.
- use the .pro file to open the project.
- Choose only the kits based on msys2/mingw
- compile.
Qt, Qt creator and msys2 are absolutly not required for compiling the project, only mingw and gcc. If someone could provide a makefile, please let me know. 
... or simply use the provided executables.

An "experimental" makefile is provided.

Dependencies:
* SDL2 : use msys2's pacman to install it.
* mpfr and GMP: They come with mingw in principle.

TODO:
A serious readme.
A serious compilation guide.
Correct spiling mistakes. :D
A lot of other things.