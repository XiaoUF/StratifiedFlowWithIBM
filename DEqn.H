
forAll(mesh.C(),i)
{
  scalar xc0;
  scalar l2;
  scalar zz;
  scalar vs;
  vector coords = mesh.C()[i];
  xc0 = xc.value() + U0.value()*runTime.value();
  l2  = (coords[0]-xc0)*(coords[0]-xc0)/Li.value()/Li.value();
  zz  = Hi.value()*Foam::exp(-l2);
  if(abs(coords[0]-xc0)>L0i.value()+0.02) {
    zz = 0.0;
  }
  D[i]  = 0.5*K.value()*(1.0-Foam::tanh(500.*(coords[2]-zz)));
  DS[i] = 0.5*KS.value()*(1.0-Foam::tanh(500.*(coords[2]-zz)));
// Calcluate the velocity
  vs = 0.5*U0.value()*(1.0-Foam::tanh(500.*(coords[2]-zz)));
  Uf[i] = vector(vs, 0.0, 0.0);

// Sponge layer
  scalar h1 = HsL.value();
  scalar h2 = HsH.value();
  if(coords[2]>h1){
    D[i]  = DSL.value()*(coords[2]-h1)/(h2-h1);
  } 

}

