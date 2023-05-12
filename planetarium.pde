int N = 500 ; 

float x[] = new float[N] ;
float y[] = new float[N] ;
float z[] = new float[N] ;
float nx[] = new float[N] ;
float ny[] = new float[N] ;
float nz[] = new float[N] ;
float dt=0.001 ;
float time=0;
float omega=1;
int skip=10;
int finishsteps=10000;

int count=0;
int iattr=0;

//Number of Attractors
int Nsims=5;

//aizawa params
float a0=0.95;
float b0=0.7;
float c0=0.65;
float d0=3.5;
float e0=0.25;
float f0=0.1;
float scale0=370;

//lorenz params
float p1=10;
float r1=28;
float b1=8/3;
float scale1=30;

//Gumowski-Mira params
float a2=0.009;
float s2=0.05;
float m2=-0.801;
float scale2=40;

//Three-Scroll Unified params
float a3=32.48;
float b3=45.84;
float c3=1.18;
float d3=0.13;
float e3=0.57;
float f3=14.7;
float scale3=4;

//Four-Wing params
float a4=0.2;
float b4=0.01;
float c4=-0.4;
float scale4=300;



//
int waittime=0;


//
// マイクの音を扱うためのライブラリをインポート
import processing.sound.*;

// マイク
AudioIn in;
// 音量を取得してくれるもの
Amplitude amp;
//



void setup(){
  fullScreen();
  //size(1000,1000) ;
  frameRate(120);
  //noSmooth() ; 
  background(0) ;
  for (int i = 0 ; i < N ; i++ ){
    float r = random(0,100);
    float theta=random(PI/100,PI/100+0.1);
    float phi=random(0,2*PI);
    x[i] = r*sin(theta)*cos(phi) ;
    y[i] = r*sin(theta)*sin(phi) ;
    z[i] = r*cos(theta) ;
  }
  // マイクを初期化（おまじない。現時点で深い意味は考えなくてOK）
  in = new AudioIn(this);
  in.start();

  // 音量の取得を開始
  amp = new Amplitude(this);
  amp.input(in);
}

void draw(){
  count+=1;
  time+=dt;
  translate(width/2,height/2);
  fill(0,0,0,50) ;
  rect(-width*2,-height*2,width*4,height*4) ; 
  fill(255);
  if (iattr==0){
    for (int step=0;step<skip;step++){
      for (int i = 0 ; i < N ; i++ ){
        nx[i] = x[i] + dt * (((z[i]+0.8*scale0)/scale0 - b0) * x[i]- d0 * y[i]);
        ny[i] = y[i] + dt * (d0 * x[i] + ((z[i]+0.8*scale0)/scale0 - b0) * y[i]);
        nz[i] = z[i] + dt * (c0 * scale0 + a0 * (z[i]+0.8*scale0) - (z[i]+0.8*scale0) * (z[i]+0.8*scale0) * (z[i]+0.8*scale0)/3/scale0/scale0 - (x[i] * x[i] + y[i] * y[i])*(1 + e0/scale0 * (z[i]+0.8*scale0))/scale0 + f0 * (z[i]+0.8*scale0) * x[i] * x[i] * x[i]/scale0/scale0/scale0 );
        stroke(100+10*randomGaussian(),249+150*randomGaussian(),255+100*randomGaussian(),100+50*randomGaussian());
        strokeWeight((2+2*atan(scale0 * y[i]*cos(omega*time)+ scale0 * x[i]*sin(omega*time))/PI));
        point((x[i]*cos(omega*time)-y[i]*sin(omega*time)),-z[i]) ;
        x[i]=nx[i];
        y[i]=ny[i];
        z[i]=nz[i];
      }
      delay(waittime);
    }
  }else if (iattr==1){
    for (int step=0;step<skip;step++){
      for (int i = 0 ; i < N ; i++ ){
        nx[i]=x[i] +p1*(y[i]-x[i])*dt/10;
        ny[i]=y[i] +(-x[i]*((z[i]+0.7)/scale1+19.3)+r1*x[i]-y[i])*dt/10;
        nz[i]=z[i] +(x[i]*y[i]/scale1-b1*((z[i]+0.7)+19.3*scale1))*dt/10;
        stroke(100+10*randomGaussian(),249+150*randomGaussian(),255+100*randomGaussian(),100+50*randomGaussian());
        strokeWeight((2+2*atan(y[i]*cos(omega*time)+x[i]*sin(omega*time))/PI));
        point((x[i]*cos(omega*time)-y[i]*sin(omega*time)),-z[i]) ;
        x[i]=nx[i];
        y[i]=ny[i];
        z[i]=nz[i]; 
      }
      delay(waittime);
    }
  }else if (iattr==2){
    for (int step=0;step<skip;step++){
      for (int i = 0 ; i < N ; i++ ){
        strokeWeight(3);
        stroke(100+10*randomGaussian(),249+100*randomGaussian(),255+100*randomGaussian(),100+50*randomGaussian());
        point((x[i]*cos(omega*time)),-(y[i]));
        nx[i]= y[i] + a2*y[i]*(1-s2*y[i]*y[i]/scale2/scale2) + m2*x[i] + 2*(1-m2)*x[i]*x[i]/scale2/(1 + x[i]*x[i]/scale2/scale2);
        ny[i]=-x[i] + m2*nx[i] + 2*(1-m2)*nx[i]*nx[i]/scale2/(1 + nx[i]*nx[i]/scale2/scale2);
        x[i]=nx[i];
        y[i]=ny[i];
        delay(waittime);
      }
      delay(waittime);
    }
  }else if (iattr==3){
    for (int step=0;step<skip;step++){
      for (int i = 0 ; i < N ; i++ ){
        x[i]=x[i] +(a3*(y[i]-x[i]) + d3*x[i]*(z[i]+430)/scale3 )*dt/10;
        y[i]=y[i] +((b3-(z[i]+430)/scale3)*x[i]+f3*y[i])*dt/10;
        z[i]=z[i] +(c3*(z[i]+430)+x[i]*y[i]/scale3-e3*x[i]*x[i]/scale3)*dt/10;
        stroke(100+10*randomGaussian(),249+150*randomGaussian(),255+100*randomGaussian(),100+50*randomGaussian());
        strokeWeight((2+2*atan(y[i]*cos(omega*time)+x[i]*sin(omega*time))/PI));
        point((x[i]*cos(omega*time)-y[i]*sin(omega*time)),-z[i]) ;
      }
      delay(waittime);
    }
  }else if (iattr==4){
    for (int step=0;step<skip;step++){
      for (int i = 0 ; i < N ; i++ ){
        x[i]=x[i] +(a4*x[i] + y[i]*z[i]/scale4)*dt*10;
        y[i]=y[i] +(b4*x[i]+c4*y[i]-x[i]*z[i]/scale4)*dt*10;
        z[i]=z[i] +(-z[i]-x[i]*y[i]/scale4)*dt*10;
        stroke(100+10*randomGaussian(),249+150*randomGaussian(),255+100*randomGaussian(),100+50*randomGaussian());
        strokeWeight((2+2*atan(y[i]*cos(omega*time)+x[i]*sin(omega*time))/PI));
        point((x[i]*cos(omega*time)-y[i]*sin(omega*time)),-z[i]) ;
      }
      delay(waittime);
    }
  }
  float A = amp.analyze();//sound import
  if (A>0.4){
    for (int i = 0 ; i < N ; i++ ){
      x[i]=x[i]+(100000*A*randomGaussian())*dt;
      y[i]=y[i]+(100000*A*randomGaussian())*dt;
      z[i]=z[i]+(100000*A*randomGaussian())*dt;
    }
  }
  if (count%finishsteps==0){
    //
    iattr=iattr+1;
    iattr=iattr%Nsims;
    count=0;
    //
    //for (int i = 0 ; i < N ; i++ ){
    //  float r = random(0,3);
    //  float theta=random(0,PI);
    //  float phi=random(0,2*PI);
    //  x[i] = r*sin(theta)*cos(phi) ;
    //  y[i] = r*sin(theta)*sin(phi) ;
    //  z[i] = r*cos(theta) ;
    //}
  }
}
