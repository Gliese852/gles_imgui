/*
 * Copyright (C) 2008  VMware, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 */

/*
 * Draw a triangle with X/EGL and OpenGL ES 2.x
 */

#define USE_FULL_GL 0

#include <iostream>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#if USE_FULL_GL
#include "gl_wrap.h"  /* use full OpenGL */
#else
#include <GLES3/gl3.h>  /* use OpenGL ES 3.x */
#endif
#include <EGL/egl.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_opengl3.h"

#include "Quaternion.h"
#include "matrix3x3.h"
#include "matrix4x4.h"
#include "vector2.h"
#include "vector3.h"
#include "Random.h"

inline constexpr double DEG2RAD(double x) { return x * (M_PI / 180.); }
inline constexpr float DEG2RAD(float x) { return x * (float(M_PI) / 180.f); }
inline constexpr double RAD2DEG(double x) { return x * (180. / M_PI); }
inline constexpr float RAD2DEG(float x) { return x * (180.f / float(M_PI)); }

struct TextureWindowData {
   uint32_t seed = 0;
   int points = 13;
   int parts = 5;
   float mindist = 3.0;
   float scale = 2.f;
   static constexpr uint32_t citySize = 384;
   uint32_t pitch = 48;
   uint32_t bsLength = citySize * pitch;

   uint8_t *bitset = nullptr;
   uint32_t *texture;
   GLuint bitsetID;

   float bars[citySize]{};
   float normalBars[citySize]{};
   float add = 0.0001f;
   int step = 10000;
   Random rand;

   uint32_t normColor = 0xFF0000FF;
   uint32_t testColor = 0xFFFF0000;
   uint32_t test1Color = 0xFF00FF00;

   float area = 0.f;

};

struct GlobalData {
   GLfloat view_rotx = 0.0;
   GLfloat view_roty = 0.0;
   GLfloat view_rotz = 0.0;
   GLfloat scalex = 0.5;
   GLfloat scaley = 0.5;
   GLfloat scalez = 0.5;
   GLint u_matrix = -1;
   GLint attr_pos = 0, attr_color = 1;

   bool show_triangle_window = true;
   bool show_texture_window = true;
   bool show_demo_window = false;

   GLuint vbo, cbo, ibo;

   TextureWindowData texture_window_data;
};

static GlobalData g;

void bitset_to_texture(uint8_t *bitset, uint32_t *texture, uint32_t citySize)
{
	std::memset(texture, 0, citySize * citySize * 4);
	for (uint32_t i = 0; i < citySize * citySize; ++i) {
		uint32_t bitset_idx = i >> 3;
		uint8_t bitmask = 1 << (i & 7);
		if (0 != (bitmask & *(bitset + bitset_idx)))
			texture[i] = 0xFFFFFFFF;
	}
}

GLuint push_texture(uint8_t *bitset, uint32_t *texture, uint32_t citySize)
{

	// bitset_to_texture(bitset, texture, citySize);

	GLuint bitsetID = 0;

	GLint last_texture;
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
	glGenTextures(1, &bitsetID);
	glBindTexture(GL_TEXTURE_2D, bitsetID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, citySize, citySize, 0, GL_RGBA, GL_UNSIGNED_BYTE, texture);
	glBindTexture(GL_TEXTURE_2D, last_texture);
	return bitsetID;
}

static void
init_texture_window(TextureWindowData &d)
{
   printf("INIT TEXTURE WINDOW\n");
   d.bitset = new uint8_t[d.bsLength];
	std::memset(d.bitset, 0, d.bsLength);
   d.texture = new uint32_t[d.citySize * d.citySize];
	std::memset(d.texture, 0, d.citySize * d.citySize * 4);
}

static void put_point(TextureWindowData &d, uint32_t x, uint32_t y, uint32_t color)
{
   uint32_t index = y * d.citySize + x;
   uint32_t &val = d.texture[index];
   val |= color;
}

static void
draw_box(TextureWindowData &d, int x0, int y0, int x1, int y1, uint32_t color)
{
   for (int x = x0; x <= x1; ++x)
      for (int y = y0; y <= y1; ++y)
         put_point(d, x, y, color);
}

// x 0..1, y 0..1
static void
draw_bar(TextureWindowData &d, int x, float y, uint32_t color)
{
   int h = y * d.citySize;
   // we will add a very rare event y = 1.0 to the last interval
   h = std::min(h, int(d.citySize - 1));
   draw_box(d, x, d.citySize - h - 1, x, d.citySize, color);
}

static void
new_value(TextureWindowData &d, float value, float array[],  uint32_t color)
{
   int x = (value + 1.f) / 2.f * d.citySize;
   array[x] += d.add;
   draw_bar(d, x, array[x], color);
   d.area += d.add;
}

static void
show_texture_window(TextureWindowData &d)
{
   ImGui::Begin("Texture window");                          // Create a window called "Hello, world!" and append into it.
   bool update = true;
   for (int i = 0; i < d.step; ++i) {
      //fixed val = d.rand.NormFixed(0, fixed(1, 1));
      fixed mean(0, 1);
      fixed maxdev(1, 1);
      fixed val = mean + maxdev * ((d.rand.Fixed() + d.rand.Fixed() + d.rand.Fixed()) * fixed(10, 15) - fixed(1, 1));
      if (val.ToFloat() > -1.f && val.ToFloat() < 1.f) new_value(d, val.ToFloat(), d.bars, d.testColor);

      //fixed val1 = d.rand.SFixed(2).Abs();
      //new_value(d, val1.ToFloat(), d.normalBars, d.test1Color);
      //double val = d.rand.Double() + d.rand.Double() + d.rand.Double();
      //val = (val / 3.0 - 0.5) * 1.9;
      //double val = d.rand.Normal(0, 0.33333);
      //if (val > -1.0 && val < 1.0) new_value(d, val, d.bars, d.testColor);
      //fixed val = d.rand.NormFixed(fixed(0, 1), fixed(1,1));
      double x = d.rand.Double_open();
      double y = d.rand.Double_open();
      double sigma = 0.33333;
      double vald = sigma * cos(2 * M_PI * x) * sqrt(-2 * log(y));

      if (vald > -1.0 && vald < 1.0) new_value(d, vald, d.normalBars, d.normColor);
      //  (* sigma (cos (* 2 pi x)) (sqrt (* -2 (log y))))
   }
   if (update) {
      d.bitsetID = push_texture(d.bitset, d.texture, d.citySize);
   }
   ImGui::SameLine();
   ImGui::Text("area: %.3f", d.area / d.citySize);
   ImTextureID my_tex_id = (ImTextureID)(intptr_t)d.bitsetID;
   //ImTextureID my_tex_id = io.Fonts->TexID;
   ImVec2 uv_min = ImVec2(0.0f, 0.0f);                 // Top-left
   ImVec2 uv_max = ImVec2(1.0f, 1.0f);                 // Lower-right
   ImVec4 tint_col = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);   // No tint
   ImVec4 border_col = ImVec4(1.0f, 1.0f, 1.0f, 0.5f); // 50% opaque white
   ImGui::Image(my_tex_id, ImVec2(d.citySize * d.scale, d.citySize * d.scale / 2), uv_min, uv_max, tint_col, border_col);
   ImGui::End();
}

static void
show_triangle_window(GlobalData &g)
{
   ImGui::Begin("Triangle control", &g.show_triangle_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
   ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
   if (ImGui::BeginTable("split", 2, ImGuiTableFlags_BordersOuter | ImGuiTableFlags_Resizable))
   {
      ImGui::AlignTextToFramePadding();

      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0); ImGui::Text("scalex");
      ImGui::TableSetColumnIndex(1); ImGui::DragFloat("##scalex", &g.scalex, 0.01f);

      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0); ImGui::Text("scaley");
      ImGui::TableSetColumnIndex(1); ImGui::DragFloat("##scaley", &g.scaley, 0.01f);

      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0); ImGui::Text("scalez");
      ImGui::TableSetColumnIndex(1); ImGui::DragFloat("##scalez", &g.scalez, 0.01f);

      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0); ImGui::Text("rotx");
      ImGui::TableSetColumnIndex(1); ImGui::DragFloat("##rotx", &g.view_rotx, 1.0f);

      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0); ImGui::Text("roty");
      ImGui::TableSetColumnIndex(1); ImGui::DragFloat("##roty", &g.view_roty, 1.0f);

      ImGui::TableNextRow();

      ImGui::TableSetColumnIndex(0); ImGui::Text("rotz");
      ImGui::TableSetColumnIndex(1); ImGui::DragFloat("##rotz", &g.view_rotz, 1.0f);

      ImGui::EndTable();
   }

   ImGui::PopStyleVar();
   ImGui::End();
}


static void
imgui_init()
{
   ImGui::CreateContext();
   ImGui_ImplOpenGL3_Init();
}


static void
imgui_prepare_data(GlobalData &g)
{
   ImGui_ImplOpenGL3_NewFrame();
   ImGui::NewFrame();
   if (g.show_demo_window) ImGui::ShowDemoWindow(&g.show_demo_window);

   show_triangle_window(g);
   show_texture_window(g.texture_window_data);
   ImGui::Render();
}


static bool
imgui_process_event(XEvent &event, Time current_time)
{
   switch (event.type) {
   case MotionNotify:
   {
      long delta = current_time - event.xmotion.time;
      if (delta > 20) return false;
      ImVec2 mouse_pos((float)event.xmotion.x, (float)event.xmotion.y);
      ImGui::GetIO().AddMousePosEvent(mouse_pos.x, mouse_pos.y);
      return true;
   }
   case ButtonPress:
   case ButtonRelease:
      int mouse_button = -1;
      switch (event.xbutton.button) {
         case 1: mouse_button = 0; break;
         case 2: mouse_button = 1; break;
         case 3: mouse_button = 2; break;

         // wheel up, down, left, right
         case 4: mouse_button = 4; break;
         case 5: mouse_button = 5; break;
         case 6: mouse_button = 6; break;
         case 7: mouse_button = 7; break;
      }
      if (mouse_button == -1)
         break;

      if (mouse_button < 4) {
         ImGui::GetIO().AddMouseButtonEvent(mouse_button, event.type == ButtonPress);
      // wheel
      } else if (event.type == ButtonPress) {
         const float speed = 0.3f;
         float wheel_y = mouse_button == 4 ? speed : mouse_button == 5 ? -speed : 0.f;
         float wheel_x = mouse_button == 6 ? speed : mouse_button == 7 ? -speed : 0.f;
         ImGui::GetIO().AddMouseWheelEvent(wheel_x, wheel_y);
         long delta = current_time - event.xmotion.time;
         if (delta > 50) return false;
      }
      return true;
   }
   return false;
}


static void
imgui_render_data()
{
   ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}


static void
imgui_destroy()
{
   ImGui::DestroyContext();
}

static void
draw(void)
{
   /* Set modelview/projection matrix */

   auto rot = matrix4x4f::RotateZMatrix(DEG2RAD(g.view_rotz));
   rot.RotateY(DEG2RAD(g.view_roty));
   rot.RotateX(DEG2RAD(g.view_rotx));

   GLint vp[4];
   glGetIntegerv(GL_VIEWPORT, vp);
   GLint width = vp[2];
   GLint height = vp[3];
   float aspect = (float)width / (float)height;

   auto scale = matrix4x4f::ScaleMatrix(g.scalex, g.scaley, g.scalez);
   auto view_scale = matrix4x4f::ScaleMatrix(1.0, aspect, 1.0);
   auto mat = view_scale * rot * scale;

   glUniformMatrix4fv(g.u_matrix, 1, GL_FALSE, mat.Data());

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   {

   //   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g.ibo);
   //   glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
   //   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

      imgui_render_data();

   }
}


/* new window size or exposure */
static void
reshape(int width, int height)
{
   glViewport(0, 0, (GLint) width, (GLint) height);
   ImGui::GetIO().DisplaySize = ImVec2((float)width, (float)height);
}


static void
make_buffers()
{
   GLfloat verts[4][2] = {
      { -1, -1 },
      {  1, -1 },
      {  1,  1 },
      { -1,  1 }
   };
   glGenBuffers(1, &g.vbo);
   glBindBuffer(GL_ARRAY_BUFFER, g.vbo);
   glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
   glEnableVertexAttribArray(g.attr_pos);
   glVertexAttribPointer(g.attr_pos, 2, GL_FLOAT, GL_FALSE, 0, 0);

   GLfloat colors[4][3] = {
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 1 },
      { 0, 1, 1 }
   };
   glGenBuffers(1, &g.cbo);
   glBindBuffer(GL_ARRAY_BUFFER, g.cbo);
   glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW);
   glEnableVertexAttribArray(g.attr_color);
   glVertexAttribPointer(g.attr_color, 3, GL_FLOAT, GL_FALSE, 0, 0);

   GLuint elements[6] = {
      0, 1, 2, 2, 3, 0
   };
   glGenBuffers(1, &g.ibo);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g.ibo);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

static void
create_shaders(void)
{
   static const char *fragShaderText =
      "precision mediump float;\n"
      "varying vec4 v_color;\n"
      "void main() {\n"
      "   gl_FragColor = v_color;\n"
      "}\n";
   static const char *vertShaderText =
      "uniform mat4 modelviewProjection;\n"
      "attribute vec4 pos;\n"
      "attribute vec4 color;\n"
      "varying vec4 v_color;\n"
      "void main() {\n"
      "   gl_Position = modelviewProjection * pos;\n"
      "   v_color = color;\n"
      "}\n";

   GLuint fragShader, vertShader, program;
   GLint stat;

   fragShader = glCreateShader(GL_FRAGMENT_SHADER);
   glShaderSource(fragShader, 1, (const char **) &fragShaderText, NULL);
   glCompileShader(fragShader);
   glGetShaderiv(fragShader, GL_COMPILE_STATUS, &stat);
   if (!stat) {
      printf("Error: fragment shader did not compile!\n");
      exit(1);
   }

   vertShader = glCreateShader(GL_VERTEX_SHADER);
   glShaderSource(vertShader, 1, (const char **) &vertShaderText, NULL);
   glCompileShader(vertShader);
   glGetShaderiv(vertShader, GL_COMPILE_STATUS, &stat);
   if (!stat) {
      printf("Error: vertex shader did not compile!\n");
      exit(1);
   }

   program = glCreateProgram();
   glAttachShader(program, fragShader);
   glAttachShader(program, vertShader);
   glLinkProgram(program);

   glGetProgramiv(program, GL_LINK_STATUS, &stat);
   if (!stat) {
      char log[1000];
      GLsizei len;
      glGetProgramInfoLog(program, 1000, &len, log);
      printf("Error: linking:\n%s\n", log);
      exit(1);
   }

   glUseProgram(program);

   if (0) {
      /* test setting attrib locations */
      glBindAttribLocation(program, g.attr_pos, "pos");
      glBindAttribLocation(program, g.attr_color, "color");
      glLinkProgram(program);  /* needed to put attribs into effect */
   }
   else {
      /* test automatic attrib locations */
      g.attr_pos = glGetAttribLocation(program, "pos");
      g.attr_color = glGetAttribLocation(program, "color");
   }

   g.u_matrix = glGetUniformLocation(program, "modelviewProjection");
   printf("Uniform modelviewProjection at %d\n", g.u_matrix);
   printf("Attrib pos at %d\n", g.attr_pos);
   printf("Attrib color at %d\n", g.attr_color);
}


static void
init(void)
{
   typedef void (*proc)();

#if 1 /* test code */
   proc p = eglGetProcAddress("glMapBufferOES");
   assert(p);
#endif

   glClearColor(0.4, 0.4, 0.4, 0.0);

   create_shaders();
   make_buffers();

}


/*
 * Create an RGB, double-buffered X window.
 * Return the window and context handles.
 */
static void
make_x_window(Display *x_dpy, EGLDisplay egl_dpy,
              const char *name,
              int x, int y, int width, int height,
              Window *winRet,
              EGLContext *ctxRet,
              EGLSurface *surfRet)
{
   static const EGLint attribs[] = {
      EGL_RED_SIZE, 1,
      EGL_GREEN_SIZE, 1,
      EGL_BLUE_SIZE, 1,
      EGL_DEPTH_SIZE, 1,
      EGL_RENDERABLE_TYPE, EGL_OPENGL_ES2_BIT,
      EGL_NONE
   };
#if USE_FULL_GL
   static const EGLint ctx_attribs[] = {
       EGL_NONE
   };
#else
   static const EGLint ctx_attribs[] = {
      EGL_CONTEXT_CLIENT_VERSION, 2,
      EGL_NONE
   };
#endif

   int scrnum;
   XSetWindowAttributes attr;
   unsigned long mask;
   Window root;
   Window win;
   XVisualInfo *visInfo, visTemplate;
   int num_visuals;
   EGLContext ctx;
   EGLConfig config;
   EGLint num_configs;
   EGLint vid;

   scrnum = DefaultScreen( x_dpy );
   root = RootWindow( x_dpy, scrnum );

   if (!eglChooseConfig( egl_dpy, attribs, &config, 1, &num_configs)) {
      printf("Error: couldn't get an EGL visual config\n");
      exit(1);
   }

   assert(config);
   assert(num_configs > 0);

   if (!eglGetConfigAttrib(egl_dpy, config, EGL_NATIVE_VISUAL_ID, &vid)) {
      printf("Error: eglGetConfigAttrib() failed\n");
      exit(1);
   }

   /* The X window visual must match the EGL config */
   visTemplate.visualid = vid;
   visInfo = XGetVisualInfo(x_dpy, VisualIDMask, &visTemplate, &num_visuals);
   if (!visInfo) {
      printf("Error: couldn't get X visual\n");
      exit(1);
   }

   /* window attributes */
   attr.background_pixel = 0;
   attr.border_pixel = 0;
   attr.colormap = XCreateColormap( x_dpy, root, visInfo->visual, AllocNone);
   attr.event_mask =
      StructureNotifyMask | ExposureMask | KeyPressMask |
      ButtonPressMask | ButtonReleaseMask | PointerMotionMask | PropertyChangeMask;
   mask = CWBackPixel | CWBorderPixel | CWColormap | CWEventMask;

   win = XCreateWindow( x_dpy, root, 0, 0, width, height,
		        0, visInfo->depth, InputOutput,
		        visInfo->visual, mask, &attr );

   /* set hints and properties */
   {
      XSizeHints sizehints;
      sizehints.x = x;
      sizehints.y = y;
      sizehints.width  = width;
      sizehints.height = height;
      sizehints.flags = USSize | USPosition;
      XSetNormalHints(x_dpy, win, &sizehints);
      XSetStandardProperties(x_dpy, win, name, name,
                              None, (char **)NULL, 0, &sizehints);
   }

#if USE_FULL_GL /* XXX fix this when eglBindAPI() works */
   eglBindAPI(EGL_OPENGL_API);
#else
   eglBindAPI(EGL_OPENGL_ES_API);
#endif

   ctx = eglCreateContext(egl_dpy, config, EGL_NO_CONTEXT, ctx_attribs );
   if (!ctx) {
      printf("Error: eglCreateContext failed\n");
      exit(1);
   }

#if !USE_FULL_GL
   /* test eglQueryContext() */
   {
      EGLint val;
      eglQueryContext(egl_dpy, ctx, EGL_CONTEXT_CLIENT_VERSION, &val);
      assert(val == 2);
   }
#endif

   *surfRet = eglCreateWindowSurface(egl_dpy, config, win, NULL);
   if (!*surfRet) {
      printf("Error: eglCreateWindowSurface failed\n");
      exit(1);
   }

   /* sanity checks */
   {
      EGLint val;
      eglQuerySurface(egl_dpy, *surfRet, EGL_WIDTH, &val);
      assert(val == width);
      eglQuerySurface(egl_dpy, *surfRet, EGL_HEIGHT, &val);
      assert(val == height);
      assert(eglGetConfigAttrib(egl_dpy, config, EGL_SURFACE_TYPE, &val));
      assert(val & EGL_WINDOW_BIT);
   }

   XFree(visInfo);

   *winRet = win;
   *ctxRet = ctx;
}

static Bool
timestamp_predicate(Display *display, XEvent  *xevent, XPointer arg)
{
   Window win = (Window)arg;
   if (xevent->type == PropertyNotify &&
         xevent->xproperty.window == win &&
         xevent->xproperty.atom == XInternAtom (display, "MY_TIMESTAMP_PROP", False)) {
      return True;
   }

   return False;
}

static Time
get_server_time (Display *dpy, Window win)
{
  const unsigned char c = 'a';
  XEvent xevent;

  Atom timestamp_prop_atom = XInternAtom(dpy, "MY_TIMESTAMP_PROP", False);

  XChangeProperty(dpy, win, timestamp_prop_atom,
        timestamp_prop_atom,
        8, PropModeReplace, &c, 1);

  XIfEvent(dpy, &xevent, timestamp_predicate, (XPointer)win);

  return xevent.xproperty.time;
}

static void
event_loop(Display *dpy, Window win, EGLDisplay egl_dpy, EGLSurface egl_surf)
{
   int redraw = 1;
   while (1) {
      XEvent event;

      XNextEvent(dpy, &event);

      switch (event.type) {
      case Expose:
         redraw = 1;
         break;
      case ConfigureNotify:
         reshape(event.xconfigure.width, event.xconfigure.height);
         break;
      case KeyPress:
         {
            char buffer[10];
            int r, code;
            code = XLookupKeysym(&event.xkey, 0);
            if (code == XK_Left) {
               g.view_roty += 5.0;
            }
            else if (code == XK_Right) {
               g.view_roty -= 5.0;
            }
            else if (code == XK_Up) {
               g.view_rotz += 5.0;
            }
            else if (code == XK_Down) {
               g.view_rotz -= 5.0;
            }
            else {
               r = XLookupString(&event.xkey, buffer, sizeof(buffer),
                                 NULL, NULL);
               if (buffer[0] == 27) {
                  /* escape */
                  return;
               }
            }
         }
         redraw = 1;
         break;
      case DestroyNotify:
         return;
      case UnmapNotify:
         break;
      default:
         if (imgui_process_event(event, get_server_time(dpy, win))) {
            redraw = 1;
         }
      }

      if (redraw) {
         imgui_prepare_data(g);
         draw();
         eglSwapBuffers(egl_dpy, egl_surf);
         redraw = 0;
      }
   }
}


static void
usage(void)
{
   printf("Usage:\n");
   printf("  -display <displayname>  set the display to run on\n");
   printf("  -info                   display OpenGL renderer info\n");
}


int
main(int argc, char *argv[])
{
   const int winWidth = 1024, winHeight = 768;
   Display *x_dpy;
   Window win;
   EGLSurface egl_surf;
   EGLContext egl_ctx;
   EGLDisplay egl_dpy;
   char *dpyName = NULL;
   GLboolean printInfo = GL_FALSE;
   EGLint egl_major, egl_minor;
   int i;
   const char *s;

   for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-display") == 0) {
         dpyName = argv[i+1];
         i++;
      }
      else if (strcmp(argv[i], "-info") == 0) {
         printInfo = GL_TRUE;
      }
      else {
         usage();
         return -1;
      }
   }

   x_dpy = XOpenDisplay(dpyName);
   if (!x_dpy) {
      printf("Error: couldn't open display %s\n",
             dpyName ? dpyName : getenv("DISPLAY"));
      return -1;
   }

   egl_dpy = eglGetDisplay(x_dpy);
   if (!egl_dpy) {
      printf("Error: eglGetDisplay() failed\n");
      return -1;
   }

   if (!eglInitialize(egl_dpy, &egl_major, &egl_minor)) {
      printf("Error: eglInitialize() failed\n");
      return -1;
   }

   s = eglQueryString(egl_dpy, EGL_VERSION);
   printf("EGL_VERSION = %s\n", s);

   s = eglQueryString(egl_dpy, EGL_VENDOR);
   printf("EGL_VENDOR = %s\n", s);

   s = eglQueryString(egl_dpy, EGL_EXTENSIONS);
   printf("EGL_EXTENSIONS = %s\n", s);

   s = eglQueryString(egl_dpy, EGL_CLIENT_APIS);
   printf("EGL_CLIENT_APIS = %s\n", s);

   make_x_window(x_dpy, egl_dpy,
                 "OpenGL ES 2.x tri", 0, 0, winWidth, winHeight,
                 &win, &egl_ctx, &egl_surf);

   XMapWindow(x_dpy, win);
   if (!eglMakeCurrent(egl_dpy, egl_surf, egl_surf, egl_ctx)) {
      printf("Error: eglMakeCurrent() failed\n");
      return -1;
   }

   if (printInfo) {
      printf("GL_RENDERER   = %s\n", (char *) glGetString(GL_RENDERER));
      printf("GL_VERSION    = %s\n", (char *) glGetString(GL_VERSION));
      printf("GL_VENDOR     = %s\n", (char *) glGetString(GL_VENDOR));
      printf("GL_EXTENSIONS = %s\n", (char *) glGetString(GL_EXTENSIONS));
   }

   imgui_init();
   init();
   init_texture_window(g.texture_window_data);

   /* Set initial projection/viewing transformation.
    * We can't be sure we'll get a ConfigureNotify event when the window
    * first appears.
    */
   reshape(winWidth, winHeight);

   event_loop(x_dpy, win, egl_dpy, egl_surf);

   eglDestroyContext(egl_dpy, egl_ctx);
   eglDestroySurface(egl_dpy, egl_surf);
   eglTerminate(egl_dpy);

   imgui_destroy();
   XDestroyWindow(x_dpy, win);
   XCloseDisplay(x_dpy);

   return 0;
}
// vim:et:sw=3:ts=3
