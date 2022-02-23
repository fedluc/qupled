

// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A TWO-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx2(int xx, int yy, int x_size) {
  return (yy * x_size) + xx;
}


// -------------------------------------------------------------------
// FUNCTION USED TO ACCESS ONE ELEMENT OF A THREE-DIMENSIONAL ARRAY
// -------------------------------------------------------------------

int idx3(int xx, int yy, int zz,
         int x_size, int y_size) {
  return (zz * x_size * y_size) + (yy * x_size) + xx;
}


// -------------------------------------------------------------------
// FUNCTION USED TO GET THE SIGN OF A NUMBER
// -------------------------------------------------------------------
int get_sign(double num) {

  if (num < 0)
    return -1;
  else
    return 1;

}  


