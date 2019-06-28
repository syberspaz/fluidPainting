final int N = 256;
final int iter = 1;
final int SCALE = 2;

int SMUDGE = 0;
float VVTHIN = 0.0000001;
float VTHIN =  0.00001;
float VTHICK = 0.0001; // bigger numbers are "thicker", smaller numbers are more fluidic

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

int IX(int x, int y){
  x = constrain(x,0,N-1);
  y = constrain(y,0,N-1);
  return x + y*N;
}

void settings() {
  size(N*SCALE,N*SCALE);
}

void setup() {
  //fluid = new Fluid(0.2,0.00000001,0.000002); // good values
  fluid = new Fluid(0.2,diffusion,viscosity);
//colorMode(HSB);
  colorRGB = new PVector(5,10,2);
  gcolorRGB = new PVector(5,10,2);
  
  img = loadImage("test.jpg");
 // fluid.setImage(img);
  cp = new ColorPicker(10,10,400,400,255);
}
void setColor(float r, float g, float b){
  colorRGB.x = r;
  colorRGB.y = g;
  colorRGB.z = b;
  gcolorRGB.set(colorRGB.x,colorRGB.y,colorRGB.z);
}

void keyPressed(){
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

void mouseReleased(){
  if(cp.isOn()) return;
  PSTRENGTH = GSTRENGTH;
      print("PSTRENGTH="+PSTRENGTH+"\n");

  //setColor(colorRGB.x,colorRGB.y,colorRGB.z);
 print(gcolorRGB);
  colorRGB = gcolorRGB;
  print(colorRGB);
  viscosity = gviscosity;
}

void mouseDragged(){
  if(cp.isOn()) return;
  int cx,cy;
  fluid.setDiffusion(diffusion);
  fluid.setViscosity(viscosity);
  cx = mouseX/SCALE;
  cy = mouseY/SCALE;
  float amtX = (mouseX - pmouseX)*0.5;
  float amtY = (mouseY - pmouseY)*0.5;
  PVector dir = new PVector(float(mouseX)-float(pmouseX),float(mouseY)-float(pmouseY));
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
    if(norm.magSq() > 0.1){
     for(int j=-BWIDTH;j<BWIDTH;j++){ 
       float xs = xx+norm.x*stepW*float(j);
       float ys = yy+norm.y*stepW*float(j);
       if(SMUDGE==0)fluid.addDensity(int(xs),int(ys),rgb);
       float dx,dy;
       dx = dir.x;//cx - xs;
       dy = dir.y;//cy - ys;
       fluid.addVelocity(int(xs),int(ys),dx*0.01,dy*0.01);
     }
    }
    //fluid.addVelocity(mouseX/SCALE,mouseY/SCALE, amtX/SCALE,amtY/SCALE);
  }
  PSTRENGTH = PSTRENGTH*0.99;
  PSTRENGTH = constrain(PSTRENGTH, 0,255);
  if(SMUDGE==0) viscosity = constrain(viscosity*1.001,VVTHIN,VTHICK);
  
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

void draw() {
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
