
class Pigment{
  float K,S;
  PVector Kv, Sv;
  
  float R,T;
  PVector Rv,Tv;
  
  float concentration;
  
  String name;
  Pigment(Pigment p){
    this.K = p.K;
    this.S = p.S;
    this.Kv = new PVector(p.Kv.x,p.Kv.y,p.Kv.z);
    this.Sv = new PVector(p.Sv.x,p.Sv.y,p.Sv.z);
    this.R = p.R;
    this.T = p.T;
    this.Rv = new PVector(p.Rv.x,p.Rv.y,p.Rv.z);
    this.Tv = new PVector(p.Tv.x,p.Tv.y,p.Tv.z);
    this.concentration = p.concentration;
    this.name = new String(p.name);
  }
  
  Pigment(){
     Kv = new PVector();
     Sv = new PVector();
     Rv = new PVector();
     Tv = new PVector();
     concentration = 0.1;
  }
  
  Pigment(String name, PVector Kv, PVector Sv, float conc){
    this.name = new String(name);
    this.Kv = new PVector(Kv.x,Kv.y,Kv.z);
    this.Sv = new PVector(Sv.x,Sv.y,Sv.z);
    this.Rv = new PVector();
    this.Tv = new PVector();
    this.concentration = conc;
  }
  
  void setName(String s){
    name = new String(s);
  }
  void setConcentration(float v) { concentration = v;}
  float getConcentration() { return concentration;}

}


final float D = 0.001;
class PAINTS{
  ArrayList<Pigment> thePaints;
  
  PAINTS(){
    thePaints = new ArrayList<Pigment>();
    thePaints.add(new Pigment("CadmiumYellow", K_CadmiumYellow, S_CadmiumYellow, D));
    thePaints.add(new Pigment("QuinacridoneRose", K_QuinacridoneRose, S_QuinacridoneRose, D));
    thePaints.add(new Pigment("FrenchUltramarine", K_FrenchUltramarine, S_FrenchUltramarine, D));
    thePaints.add(new Pigment("CeruleanBlue", K_CeruleanBlue, S_CeruleanBlue, D));
    thePaints.add(new Pigment("HookersGreen", K_HookersGreen, S_HookersGreen, D));
    thePaints.add(new Pigment("HansaYellow", K_HansaYellow, S_HansaYellow, D));
    thePaints.add(new Pigment("BrilliantOrange", K_BrilliantOrange, S_BrilliantOrange, D));
    thePaints.add(new Pigment("CadmiumRed", K_CadmiumRed, S_CadmiumRed, D));
    thePaints.add(new Pigment("IndianRed", K_IndianRed, S_IndianRed, D));
    thePaints.add(new Pigment("InterferenceLilac", K_InterferenceLilac, S_InterferenceLilac, D));
    thePaints.add(new Pigment("PhthaloGreen", K_PhthaloGreen, S_PhthaloGreen, D));
    thePaints.add(new Pigment("BurntUmber", K_BurntUmber, S_BurntUmber, D));
  } 
  Pigment get(String name){
    for(int i=0;i<thePaints.size();i++){
      Pigment p = thePaints.get(i);
      if(p.name.equals(name)) return p;
    }
    return null;
  }

}



// this is a pixel that contains paint
// each pixel will contain a list of pigments
class PaintPixel{
  ArrayList<Pigment> pigments;

  float volume; // current amount of paint in this pixel
  
  PaintPixel(){
    pigments = new ArrayList<Pigment>();
    volume = 0.0;
  }
  
  void advectOut(float vx, float vy, float dt){
    float outX = constrain(1.0f - Math.abs(vx*dt),0.01,20000.f);
    float outY = constrain(1.0f - Math.abs(vy*dt),0.01,20000.f);
    float totalOut = constrain(outX*outY,0.01,1.0);
    float vol = volume*totalOut;
    for(int i=0;i<pigments.size();i++){
      Pigment p = pigments.get(i);
      p.concentration *= totalOut;
      
    }
    volume = vol;
  }
  void addPaint(Pigment p){
     int ind;
     if((ind = findPigment(p.name))>=0){
        /* found the pigment within this list */
        Pigment pp = pigments.get(ind);
        pp.concentration += p.concentration;
        this.volume = pp.concentration;
      }
      else {
        Pigment pp = new Pigment(p);
        pigments.add(pp);
        this.volume += pp.concentration;
      }
  }
  void advectIn(float vx, float vy, float dt, PaintPixel pX, PaintPixel pY){
    advectInPixel(vx*dt,pX);
    advectInPixel(vy*dt,pY);
  }
  
  int findPigment(String name){
    for(int i=0;i<pigments.size();i++){
      if(pigments.get(i).name.equals(name)) return i;
    }
    return -1;
  }
  void advectInPixel(float delta, PaintPixel pix){
    for(int i=0;i<pix.pigments.size();i++){
      Pigment p = pix.pigments.get(i);
      int ind=0;
      if((ind = findPigment(p.name))>=0){
        /* found the pigment within this list */
        Pigment pp = pigments.get(ind);
        pp.concentration += delta;
      }
      else {
        Pigment pp = new Pigment(p);
        pigments.add(p);
        pp.concentration +=delta;  // the original pigment concentration from incoming pixel has already been copied in the copy constructor
      }
      volume += pix.volume*delta;
    }
  
  }

}
