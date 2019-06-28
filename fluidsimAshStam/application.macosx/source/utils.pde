

PVector rybToRgb(float r, float y, float b){
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



PVector rgbToRyb (float r,float g,float b) {
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
