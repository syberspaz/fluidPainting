import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class fluidsimAshStam extends PApplet {

final int N = 256;
final int iter = 1;
final int SCALE = 2;

int SMUDGE = 0;
float VVTHIN = 0.0000001f;
float VTHIN =  0.00001f;
float VTHICK = 0.0001f; // bigger numbers are "thicker", smaller numbers are more fluidic

float diffusion = 0;//0.00000001;
float viscosity = VTHIN; // bigger numbers are "thicker", smaller numbers are more fluidic
float gviscosity = VTHIN;

Fluid fluid;
float t = 0;
int brushSize = 7; // 2 = 4x4
PVector colorRGB,gcolorRGB;
float r=0;
float g = 0;
float b = 0;

PImage img;
ColorPicker cp;

float PSTRENGTH=20;
float GSTRENGTH=20;

public int IX(int x, int y){
  x = constrain(x,0,N-1);
  y = constrain(y,0,N-1);
  return x + y*N;
}

public void settings() {
  size(N*SCALE,N*SCALE);
}

public void setup() {
  //fluid = new Fluid(0.2,0.00000001,0.000002); // good values
  fluid = new Fluid(0.2f,diffusion,viscosity);
//colorMode(HSB);
  colorRGB = new PVector(5,10,2);
  gcolorRGB = new PVector(5,10,2);
  
  img = loadImage("test.jpg");
 // fluid.setImage(img);
  cp = new ColorPicker(10,10,400,400,255);
}
public void setColor(float r, float g, float b){
  colorRGB.x = r;
  colorRGB.y = g;
  colorRGB.z = b;
  gcolorRGB.set(colorRGB.x,colorRGB.y,colorRGB.z);
}

public void keyPressed(){
  if(key == '1'){
    setColor(1,0,0);
     SMUDGE = 0;
  }
  else if(key == '2'){
    setColor(0,1,0);
     SMUDGE = 0;//SMUDGE==0?1:0;
  }
  else if(key == '3'){
    setColor(0,0,1);
     SMUDGE = 0;//SMUDGE==0?1:0;
  }
  else if(key == '3'){
    setColor(0,1,1);
     SMUDGE = 0;//SMUDGE==0?1:0;
  }
  else if(key == 'c'){
    cp.toggle();
  // fluid.addVelocityField();
  }
  else if(key == 'd'){
   diffusion /=2.f;
      print("diffusion="+diffusion+"\n");
  }
  else if(key == 'D'){
   diffusion *=2.f;
   print("diffusion="+diffusion+"\n");

  }
  else if(key == 'v'){
   viscosity /=2.f;
      print("viscosity="+viscosity+"\n");

  }
  else if(key=='V'){
   viscosity *=2.f;
   print("viscosity="+viscosity+"\n");
  }
  else if(key == 's'){
    setColor(1,1,1);
  SMUDGE = SMUDGE==0?1:0;
  viscosity = VTHIN;
  }
  else if(key == 'b'){
    brushSize--;
    brushSize = constrain(brushSize,1,20);
  }
  else if(key == 'B'){
    brushSize++;
    brushSize = constrain(brushSize,1,20);

  }

  else if (key == 'p'){
     GSTRENGTH = PSTRENGTH--;
    print("PSTRENGTH="+PSTRENGTH+"\n");
  }
  else if(key == 'P'){
  GSTRENGTH = PSTRENGTH++;
    print("PSTRENGTH="+PSTRENGTH+"\n");

}
}

public void mouseReleased(){
  if(cp.isOn()) return;
  PSTRENGTH = GSTRENGTH;
      print("PSTRENGTH="+PSTRENGTH+"\n");

  //setColor(colorRGB.x,colorRGB.y,colorRGB.z);
 print(gcolorRGB);
  colorRGB = gcolorRGB;
  print(colorRGB);
  viscosity = gviscosity;
}

public void mouseDragged(){
  if(cp.isOn()) return;
  int cx,cy;
  fluid.setDiffusion(diffusion);
  fluid.setViscosity(viscosity);
  cx = mouseX/SCALE;
  cy = mouseY/SCALE;
  float amtX = (mouseX - pmouseX)*0.5f;
  float amtY = (mouseY - pmouseY)*0.5f;
  PVector dir = new PVector(PApplet.parseFloat(mouseX)-PApplet.parseFloat(pmouseX),PApplet.parseFloat(mouseY)-PApplet.parseFloat(pmouseY));
  dir.normalize();

  int BLENGTH=brushSize;
  int BWIDTH=brushSize;
  float stepL = BLENGTH/5;
  float stepW = BWIDTH/5;
 // float s = NUMS/(brushSize);
  PVector norm = new PVector(-dir.y, dir.x);
  PVector rgb = new PVector();//.normalize();
  rgb.set(colorRGB.x*PSTRENGTH, colorRGB.y*PSTRENGTH, colorRGB.z*PSTRENGTH);
  for(int i=0;i<BLENGTH;i++){
    float xx = cx - dir.x*stepL*i;
    float yy = cy - dir.y*stepL*i;
    if(norm.magSq() > 0.1f){
     for(int j=-BWIDTH;j<BWIDTH;j++){ 
       float xs = xx+norm.x*stepW*PApplet.parseFloat(j);
       float ys = yy+norm.y*stepW*PApplet.parseFloat(j);
       if(SMUDGE==0)fluid.addDensity(PApplet.parseInt(xs),PApplet.parseInt(ys),rgb);
       float dx,dy;
       dx = dir.x;//cx - xs;
       dy = dir.y;//cy - ys;
       fluid.addVelocity(PApplet.parseInt(xs),PApplet.parseInt(ys),dx*0.01f,dy*0.01f);
     }
    }
    //fluid.addVelocity(mouseX/SCALE,mouseY/SCALE, amtX/SCALE,amtY/SCALE);
  }
  PSTRENGTH = PSTRENGTH*0.99f;
  PSTRENGTH = constrain(PSTRENGTH, 0,255);
  if(SMUDGE==0) viscosity = constrain(viscosity*1.001f,VVTHIN,VTHICK);
  
}
  
//  for(int i=-brushSize;i<=brushSize;i++)
//  for(int j=-brushSize;j<=brushSize;j++)
//  {
//    if(SMUDGE==1) {
      
//      fluid.addVelocity(cx,cy, amtX*0.0001,amtY*0.0001);
//    }
//    else fluid.addDensity(cx+i,cy+j,colorRGB);
//  }
//  fluid.addVelocity(mouseX/SCALE,mouseY/SCALE, amtX/SCALE,amtY/SCALE);
//}

public void draw() {
background(0);
  if(cp.isOn())
  {
  cp.render();
 int c = cp.getColor();
  float r = (red(c))/255.f;
  float g = (green(c))/255.f;
  float b = (blue(c))/255.f;

  setColor(r,g,b);
  }else
  {
      
  fluid.renderD();
  fluid.step();
  noFill();
  stroke(color(gcolorRGB.x*100,gcolorRGB.y*100,gcolorRGB.z*100));
  strokeWeight(brushSize);
  square(mouseX-brushSize*2, mouseY-brushSize*2, brushSize*4);
  }
}


public class ColorPicker 
{
  int x, y, w, h, c;
  PImage cpImage;
  boolean isON = false;
  public ColorPicker ( int x, int y, int w, int h, int c )
  {
    this.x = x;
    this.y = y;
    this.w = w;
    this.h = h;
    this.c = c;
    
    cpImage = new PImage( w, h );
    
    init();
  }
  public void toggle(){
    isON = (isON==true)?false:true;
  }
  public boolean isOn(){
    return isON;
  }
  
  public int getColor(){
  return this.c;
  }
  
  private void init ()
  {
    // draw color.
    int cw = w - 60;
    for( int i=0; i<cw; i++ ) 
    {
      float nColorPercent = i / (float)cw;
      float rad = (-360 * nColorPercent) * (PI / 180);
      int nR = (int)(cos(rad) * 127 + 128) << 16;
      int nG = (int)(cos(rad + 2 * PI / 3) * 127 + 128) << 8;
      int nB = (int)(Math.cos(rad + 4 * PI / 3) * 127 + 128);
      int nColor = nR | nG | nB;
      
      setGradient( i, 0, 1, h/2, 0xFFFFFF, nColor );
      setGradient( i, (h/2), 1, h/2, nColor, 0x000000 );
    }
    
    // draw black/white.
    drawRect( cw, 0,   30, h/2, 0xFFFFFF );
    drawRect( cw, h/2, 30, h/2, 0 );
    
    // draw grey scale.
    for( int j=0; j<h; j++ )
    {
      int g = 255 - (int)(j/(float)(h-1) * 255 );
      drawRect( w-30, j, 30, 1, color( g, g, g ) );
    }
  }

  private void setGradient(int x, int y, float w, float h, int c1, int c2 )
  {
    float deltaR = red(c2) - red(c1);
    float deltaG = green(c2) - green(c1);
    float deltaB = blue(c2) - blue(c1);

    for (int j = y; j<(y+h); j++)
    {
      int c = color( red(c1)+(j-y)*(deltaR/h), green(c1)+(j-y)*(deltaG/h), blue(c1)+(j-y)*(deltaB/h) );
      cpImage.set( x, j, c );
    }
  }
  
  private void drawRect( int rx, int ry, int rw, int rh, int rc )
  {
    for(int i=rx; i<rx+rw; i++) 
    {
      for(int j=ry; j<ry+rh; j++) 
      {
        cpImage.set( i, j, rc );
      }
    }
  }
  
  public void render ()
  {
    image( cpImage, x, y );
    if( mousePressed &&
  mouseX >= x && 
  mouseX < x + w &&
  mouseY >= y &&
  mouseY < y + h )
    {
      c = get( mouseX, mouseY );
    }
    fill( c );
    rect( x, y+h+10, 20, 20 );
  }
  
}

float tt=0;
class Fluid {
    int size;
    float dt;
    float diff;
    float visc;
    
    float[] sR;
    float[] sG;
    float[] sB;
    float[] densityR;
    float[] densityG;  // testing: this will actually be Y, using RYB now
    float[] densityB;
    
    float[] Vx;
    float[] Vy;

    float[] Vx0;
    float[] Vy0;

public void renderD(){
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
       float x = i*SCALE;
       float y = j*SCALE;
       
       float r = this.densityR[IX(i,j)];
       float Y = this.densityG[IX(i,j)];
       float b = this.densityB[IX(i,j)];
       PVector rgb = rybToRgb(r,Y,b);
              
       fill(rgb.x,rgb.y,rgb.z);
       noStroke();
       square(x,y,SCALE);
    }
  } 
}
public void fadeD(){
 // for(int i=0;i<this.densityR.length;i++){
//    float d = density[i];
//    density[i] = constrain(d-0.01,0,255);
 // }
  
}


    Fluid(float dt, float diffusion, float viscosity) {
   
    this.size = N;
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;
    
    this.sR = new float[N*N];
    this.sG = new float[N*N];
    this.sB = new float[N*N];

    this.densityR = new float[N*N];
    this.densityG = new float[N*N];
    this.densityB = new float[N*N];
    
    this.Vx = new float[N*N];
    this.Vy = new float[N*N];
    
    this.Vx0 = new float[N*N];
    this.Vy0 = new float[N*N];
    
    }
    // add some weird random velocity thing
    public void addVelocityField(){
      float cx,cy;
      cx = 0.5f*width;
      cy = 0.5f*height;
      
      for(int j=0;j<N;j++){
      for(int i=0;i<N;i++)
      {
        float vx = this.Vx[IX(i,j)];
        float vy = this.Vy[IX(i,j)];        
        float amtX = vx*(cx-i);
        float amtY = vy*(cy-j);
            this.Vx[IX(i,j)] += amtX;
            this.Vy[IX(i,j)] += amtY;

      }
    }
  }
  public void setImage(PImage img){
      PVector rgb = new PVector();
 
    for(int j=0;j<N;j++){
    for(int i=0;i<N;i++){
      int c = img.get(i*SCALE,j*SCALE);
      rgb.set(red(c), green(c), blue(c));
      this.addDensity(i,j,rgb);
    }
    }
    
  }
  
    public void setViscosity(float v){
    this.visc = v;
    }
    public void setDiffusion(float d){
    this.diff = d;
    }
    
    public void addDensity(int x, int y, PVector rgb){
      PVector ryb = rgbToRyb(rgb.x,rgb.y,rgb.z);
      
      this.densityR[IX(x,y)] += ryb.x;
      this.densityG[IX(x,y)] += ryb.y;
      this.densityB[IX(x,y)] += ryb.z;
      
    }
    
    public void addVelocity(int x, int y, float amountX, float amountY) {
      int index = IX(x,y);
      this.Vx[index] += amountX;
      this.Vy[index] += amountY;
    }
    
    
    public void step(){
      float visc     = this.visc;
      float diff     = this.diff;
      float dt       = this.dt;
      float[] Vx      = this.Vx;
      float[] Vy      = this.Vy;
      float[] Vx0     = this.Vx0;
      float[] Vy0     = this.Vy0;
      float[] sR       = this.sR;
      float[] sG       = this.sG;
      float[] sB       = this.sB;
       
      float[] densityR = this.densityR;
      float[] densityG = this.densityG;
      float[] densityB = this.densityB;

      
      diffuse(1, Vx0, Vx, visc,dt);
      diffuse(2, Vy0, Vy, visc,dt);
      
      project(Vx0, Vy0, Vx, Vy);
      
      advect(1, Vx, Vx0, Vx0, Vy0, dt);
      advect(2, Vy, Vy0, Vx0, Vy0, dt);
      
      project(Vx, Vy,  Vx0, Vy0 );
      
      diffuse(0, sR, densityR, diff, dt);
      diffuse(0, sG, densityG, diff, dt);
      diffuse(0, sB, densityB, diff, dt);
      
      advect(0, densityR, sR, Vx, Vy,dt);
      advect(0, densityG, sG, Vx, Vy,dt);
      advect(0, densityB, sB, Vx, Vy,dt);
     
    }
}
public void set_bnd(int b, float[] x)
{
      for(int i = 1; i < N - 1; i++) {
          x[IX(i, 0 )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1  )];
          x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
      }
  
      for(int j = 1; j < N - 1; j++) {
          x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
          x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
      }
    x[IX(0, 0)]     = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)]   = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);
    x[IX(N-1, 0)]   = 0.5f * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);    
}

public void lin_solve(int b, float[] x, float[] x0, float a, float c)
{
    float cRecip = 1.0f / c;
    for (int k = 0; k < iter; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a*(    x[IX(i+1, j  )]
                                    +x[IX(i-1, j )]
                                    +x[IX(i  , j+1 )]
                                    +x[IX(i  , j-1 )]
                             
                           )) * cRecip;
                }
            }
        
        set_bnd(b, x);
    }
}


public void diffuse (int b, float[] x, float[] x0, float diff, float dt)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

public void project(float []velocX, float[] velocY, float[] p, float[] div)
{
 
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                         velocX[IX(i+1, j  )]
                        -velocX[IX(i-1, j   )]
                        +velocY[IX(i  , j+1 )]
                        -velocY[IX(i  , j-1  )]
                    )/N;
                p[IX(i, j)] = 0;
            }
        }
    
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);
    
  
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N;
              
            }
        }
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

public void advect(int b, float[] d, float[] d0,  float[] velocX, float[] velocY,float dt)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floor(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floor(y);
                j1 = j0 + 1.0f; 
              
               
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
               
                
                int i0i = PApplet.parseInt(i0);
                int i1i = PApplet.parseInt(i1);
                int j0i = PApplet.parseInt(j0);
                int j1i = PApplet.parseInt(j1);
              

                d[IX(i, j)] = 
                    s0 * (t0*d0[IX(i0i, j0i)]+ t1*d0[IX(i0i, j1i)]) + 
                    s1 * (t0*d0[IX(i1i, j0i)]+ t1*d0[IX(i1i, j1i)]);
            }
        }
    
    set_bnd(b, d);
}


public PVector rybToRgb(float r, float y, float b){
  float iW = min(r,y,b);
  float R = r - iW;
  float Y = y - iW;
  float B = b - iW;
  
  float mY = max(R,Y,B);
  float G = min(Y,B);
  
  Y -= G;
  B -= G;
  
  if(B > 0 && G > 0){
    B *= 2.0f;
    G *= 2.0f;
  }
  
  R += Y;
  G += Y;
  
  float mG = max(R,G,B);
  if(mG > 0){
    float N = mY / mG;
    R *= N;
    G *= N;
    B *= N;
  }
  
  R += iW;
  G += iW;
  B += iW;
  
  return new PVector(R,G,B);
}



public PVector rgbToRyb (float r,float g,float b) {
  float iW = min(r,g,b);
  float R = r - iW;
  float G = g - iW;
  float B = b - iW;
  
  float mG = max(R,G,B);
  float iY = min(R,G);
  R -=iY;
  G -= iY;
  if(B >0 && G > 0){
    B /=2.f;
    G /=2.f;
  }
  
  iY += G;
  B += G;
  
  float mY = max(R,iY,B);
  if(mY >0){
    float N = mG/mY;
    R *= N;
    iY *= N;
    B *= N;
  }
  R += iW;
  iY += iW;
  B += iW;
  return new PVector(R,iY,B);
}
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "fluidsimAshStam" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
