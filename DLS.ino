#include "HX711.h"

const int dout1 = 2;
const int sck1 = 3;
const int dout2 = 4;
const int sck2 = 5;
int inc = 0;

HX711 scale1;
HX711 scale2;

void setup() {
  Serial.begin(57600);
  scale1.begin(dout1, sck1);
  scale2.begin(dout2, sck2);
  scale1.set_scale(1027.540610);
  scale2.set_scale(-995.016598);
  scale1.tare();
  scale2.tare();
}

void loop() {
  //Serial.print(inc);
  Serial.print("  ");
  //Serial.print(millis());
  Serial.print("  ");
  //Serial.print("Lift: ");
  Serial.print(scale1.get_units(10), 5);
  //Serial.print("      Drag: ");
  Serial.print("  ");
  Serial.println(scale2.get_units(10), 5);
  inc = inc + 1;
}