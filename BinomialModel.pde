import interfascia.*;

float[][] angleGens =
  { /// fair approximately, where x is an oriented axis cycle matrix (x^(2n) = -1, or x^(2n+1) = 1), the coefficients a_k of a solution (log x) = Σ a_k x^k 
    {0}, // we could do these more precisely in code, but instead we got them from Maxima; five digits ought to be enough for anyone...
    {0, 1.0},
    {0,-1.00000, 1.00000},
    {0, 1.41421,-1.00000, 1.41421},
    {0,-0.85065, 0.52573,-0.52573, 0.85065},
    {0, 2.00000,-1.15470, 1.00000,-1.15470, 2.00000},
    {0,-1.15238, 0.63952,-0.51286, 0.51286,-0.63952, 1.15238},
    {0, 2.61312,-1.41421, 1.08239,-1.00000, 1.08239,-1.41421, 2.61312},
    {0,-1.46190, 0.77786,-0.57735, 0.50771,-0.50771, 0.57735,-0.77786, 1.46190}
  };
    
class dividedSimplex {
  
  int theDimension;
  int corners;
  int transSign;
  float subdivRatio;
  float theAxesX[]; 
  float theAxesY[]; 
  float centroidX, centroidY;
  float theGensX[];
  float theGensY[];
  float theAngles[];
  float lastVertexScale;
  float realCx,realCy;

  public
  dividedSimplex(int dimension){
    theDimension=dimension;
    transSign = ((theDimension & 1) == 0 ) ? -1 : 1 ;
    corners = 1 + dimension;
    subdivRatio = .7;
    lastVertexScale = sqrt(corners) -1;
    theAxesX = new float[corners];
    theGensX = new float[dsDim];
    theAxesY = new float[corners];
    theGensY = new float[dsDim];
    theAngles = new float[dsDim];
    centroidX = 0;
    centroidY = 0;
    for ( int j = 0 ; j < dsDim ; j++ ){
      theGensX[j] = random(0,1);
      theGensY[j] = random(0,1);
      centroidX += theGensX[j] / dsDim;
      centroidY += theGensY[j] / dsDim;
    }
    
  }
  
  void renorm(){
    float nx = 0 ; // as r;
    float viewSize = min(height,width);
    int j;
    for ( j = 0 ; j < dsDim; j++ ){
      nx += theGensX[j]*theGensX[j];
    }
    nx = sqrt(nx);
    centroidX = 0;
    for (j = 0; j < dsDim ; j++ ){
      theGensX[j] /= nx;
      centroidX += theGensX[j]/dsDim;
    }
    nx = 0; // as inner-product <axesX . axesY>
    for (j = 0; j < dsDim ; j++ ){
      nx += theGensX[j] * theGensY[j];
    }
    for (j = 0; j < dsDim ; j++ ){
      theGensY[j] -= nx * theGensX[j];
    }
    nx = 0;
    for (j = 0; j < dsDim ; j++ ){
      nx += theGensY[j] * theGensY[j];
    }
    nx = sqrt(nx);
    centroidY = 0;
    for (j = 0; j < dsDim ; j++ ){
      theGensY[j] /= nx;
      centroidY += theGensY[j]/dsDim;
    }
    
    for (j = 0 ; j < dsDim ; j++ ) {
      theAxesX[j] = ( width + viewSize * theGensX[j] )/ 2;
      theAxesY[j] = ( height + viewSize * theGensY[j] ) / 2;
    }
    theAxesX[dsDim] = ( width - viewSize * centroidX * lastVertexScale ) / 2 ; 
    theAxesY[dsDim] = ( height - viewSize * centroidY * lastVertexScale ) / 2 ;
  }

  void anglePerturb(){
  /* The angles, a short list of Euler Angles, should be Small, about two degrees.
  We perturb them by about half a degree
  */
    float resize = (turnPerFrame)/(turnPerFrame + bumpPerFrame);
    float nx, ny, s, c;
    int j = 0 ;
  
    for (j = 0 ; j < dsDim ; j++ ){
       theAngles[j] += (random(2) - 1) * bumpPerFrame ; // previously:  (rand(2)-1)/ 600 ; 
       theAngles[j] *= resize; // previously , resize = .95 , and turnPerFrame was left implicit // *=  resize ;
    }
    for ( j = 0; j < dsDim-1 ; j++){
      s = sin(theAngles[j]);
      c = cos(theAngles[j]);
      nx = theGensX[j] * c + theGensX[j+1] * s;
      ny = theGensX[j+1] * c - theGensX[j] * s;
      theGensX[j] = nx; theGensX[j+1] = ny;
      nx = theGensY[j] * c + theGensY[j+1]*s;
      ny = theGensY[j+1]*c - theGensY[j] * s;
      theGensY[j] = nx;
      theGensY[j+1] = ny;
    }
    s = sin(theAngles[dsDim-1]) ; c = cos(theAngles[dsDim-1]);
    nx = theGensX[j] * c + theGensX[0] * s;
    ny = theGensX[0] * c - theGensX[j] * s;
    theGensX[j] = nx; theGensX[0] = ny;
    nx = theGensY[j] * c + theGensY[0] * s;
    ny = theGensY[0] * c - theGensY[j] * s;
    theGensY[j] = nx; theGensY[0] = ny;
    
    this.renorm();

  }

  void alternatePerturb(){
    float angleOne = turnPerFrame * (mouseX - width/2)/width;
    float angleTwo = turnPerFrame * (mouseY - height/2)/height/theDimension;
    float newGens[] = new float[theDimension];
    
    float c1 = cos(angleOne);
    float s1 = sin(angleOne);
    
    int j,k;
    for ( j = 0 ; j < theDimension ; j++ ){
      if ( j < 2 ) newGens[j] = c1 * theGensX[j] + (1-2*j) * s1 * theGensX[1-j];
      else newGens[j] = theGensX[j];
    }
    for ( j = 0 ; j < theDimension ; j++ ){
      theGensX[j] = newGens[j];
      for ( k = j+1 ; k < theDimension ; k++ ){
        theGensX[j] += angleTwo * newGens[k] * (angleGens[theDimension-1][k-j]);
      }
      for ( k = j-1 ; k >= 0 ; k-- ){
        theGensX[j] += transSign * angleTwo * newGens[k] * angleGens[theDimension-1][k+theDimension - j];
      }
    }
    // and again, for the ys.
    for ( j = 0 ; j < theDimension ; j++ ){
      if ( j < 2 ) newGens[j] = c1 * theGensY[j] + (1-2*j) * s1 * theGensY[1-j];
      else newGens[j] = theGensY[j];
    }
    for ( j = 0 ; j < theDimension ; j++ ){
      theGensY[j] = newGens[j];
      for ( k = j+1 ; k < theDimension ; k++ ){
        theGensY[j] += angleTwo * newGens[k] * (angleGens[theDimension-1][k-j]);
      }
      for ( k = j-1 ; k >= 0 ; k-- ){
        theGensY[j] += transSign * angleTwo * newGens[k] * angleGens[theDimension-1][k + theDimension - j];
      }
    }
    this.renorm();
  }
 
  void drawMe_joinedUp(){
    // we only really want to draw the *edges* of the divided simplex;  
    // there are corners_choose_two long edges joining proper vertices, and 2 * corners_choose_three
    //  short edges joining various of the long edges "centers" ...
    // ... though we may try adjusting the ratio of the subdivisions.
    // 
    int l, m, p;
    for (l = 0 ; l < corners ; l++ ){
      for ( m = l + 1 ; m < corners ; m++){
        // draw a long edge;
        line(theAxesX[l],theAxesY[l],theAxesX[m],theAxesY[m]);
        for ( p = m + 1 ; p < corners ; p++ ){
        // draw the short edges! IN the face [l<m<p], join the midpoints of [l<m] and [m<p] to that of [l<p]
        // toDo: decide What Colors
          float x1 = subdivRatio * theAxesX[l] + (1-subdivRatio)*theAxesX[p];
          float y1 = subdivRatio * theAxesY[l] + (1-subdivRatio)*theAxesY[p];
          float x2 = subdivRatio * theAxesX[l] + (1-subdivRatio)*theAxesX[m];
          float y2 = subdivRatio * theAxesY[l] + (1-subdivRatio)*theAxesY[m]; 
          line(x1,y1,x2,y2);
          x2 = subdivRatio * theAxesX[m] + (1-subdivRatio) * theAxesX[p] ;
          y2 = subdivRatio * theAxesY[m] + (1-subdivRatio) * theAxesY[p] ;
          line(x1,y1,x2,y2);
        }
      }
    }
  }
  
  void drawMe_blocks(){
    int l, m, n, z;
    float x1 , x2 , y1 , y2;
    float s = 0.35;
    float t = 0.55;
    float eps = 0.1;
    colorMode(HSB,100);
    for(l = 0; l < corners ; l++){         // I Really Wish there were a Less-Loop-Nesting-Way of doing this...
      for (m = l ; m < corners ; m++ ){
 //       stroke( (61 * m)%100 , 100 , 85 );
        for ( n = m ; n < corners ; n++){
          for(z = l+1 ; z <= m ; z++) { // l < z <= m <= n
//          stroke ( (61 * z + 53 * l) %100 , 100 , 85 ); // opt B
            stroke ( (31 * (2*m+1) ) %100 , 100, 85 );  // opt C
            x1 = s * theAxesX[l] + t * theAxesX[n] + eps * theAxesX[m];
            y1 = s * theAxesY[l] + t * theAxesY[n] + eps * theAxesY[m];
            x2 = s * theAxesX[z] + t * theAxesX[n] + eps * theAxesX[m];
            y2 = s * theAxesY[z] + t * theAxesY[n] + eps * theAxesY[m];
            line(x1,y1,x2,y2);}
          for(z = m ; z < n ; z++) { // l <= m <= z < n
//            stroke( (61 * n + 53 * z) % 100 , 100 , 85 ); // opt B
            stroke ( (31 * 2 * m )%100 , 100, 85 ) ; // opt C 
            x1 = s * theAxesX[l] + t * theAxesX[z] + eps * theAxesX[m];
            y1 = s * theAxesY[l] + t * theAxesY[z] + eps * theAxesY[m];
            x2 = s * theAxesX[l] + t * theAxesX[n] + eps * theAxesX[m];
            y2 = s * theAxesY[l] + t * theAxesY[n] + eps * theAxesY[m];
            line(x1,y1,x2,y2);
          }
        }
      }
    }
    colorMode(RGB,255);
    stroke(0);
  }
}

float turnPerFrame = 3.0 / 57.3; // about three degrees  
float bumpPerFrame = 0.10 / 300.0; // about 1/60th of a degree... 


dividedSimplex ds[];
int dsDim;
int corners;
IFButton buttonDecr,buttonIncr;
GUIController controller;

void setup(){
  size(640,480);
  // smooth(4);
  noFill();
  
  surface.setTitle("Binomial Illustrated");
  surface.setResizable(true);
  
  dsDim = 4;
  
  controller = new GUIController(this);
  
  buttonDecr = new IFButton("<",0,0);
  buttonIncr = new IFButton(">", buttonDecr.getWidth(),0);
  controller.add(buttonDecr);
  controller.add(buttonIncr);
  
  buttonDecr.addActionListener(this);
  buttonIncr.addActionListener(this);
  
  corners = dsDim + 1 ;
  
  ds = new dividedSimplex[10];
  
  ds[dsDim] = new dividedSimplex(dsDim);
}

void draw(){
  background(128);
  if ( (mouseX | mouseY) == 0 ) ds[dsDim].anglePerturb();
  else ds[dsDim].alternatePerturb();
  /* for ( int j = 0 ; j < 5 ; j++ ){
    theAxesX[j] = (random(1) < .5) ? random(theAxesX[j],theAxesX[j]*.99 + height * .01 ) : random(theAxesX[j]*.99, theAxesX[j]);
    theAxesY[j] = (random(1) < .5) ? random(theAxesY[j],theAxesY[j]*.99 + height * .01 ) : random(theAxesY[j]*.99, theAxesY[j]);
  } */
  ds[dsDim].drawMe_blocks();
  // saveFrame("binomialIllustrn_####.png");
}

void actionPerformed(GUIEvent e){
  Object src = e.getSource();
  if ( src == buttonDecr){
    if ( dsDim > 2 ) dsDim-- ;
    if (null == ds[dsDim]) ds[dsDim] = new dividedSimplex(dsDim) ;
  }
  else if (src == buttonIncr){
    if ( dsDim < 9 ) dsDim++;
    if (null == ds[dsDim]) ds[dsDim] = new dividedSimplex(dsDim);
  }
}
