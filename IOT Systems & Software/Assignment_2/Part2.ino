/*
#Ref
https://www.makerguides.com/interfacing-esp32-and-tcs34725-rgb-color-sensor/
https://github.com/adafruit/Adafruit_TCS34725

#Wiring
#RGB sensor
TCS34725 - ESP32
3v3      - 3v3
SCL      - P22  
SDA      - P21
GND      - GND

#On board LED
ESP32 - PIN2

GREEN LED - 25
RED LED - 33
YELLOW LED - 32
*/

#include "Wire.h"
#include "Adafruit_TCS34725.h"

#define MAX_THRESHOLD 100
#define MIN_THRESHOLD 90
#define YELLOW_LED 32
#define RED_LED 33
#define GREEN_LED 25
// #define ONBOARD_LED 2
#define TWO_SECONDS 2000 // 2000 ms = 2 second
#define ONE_SECOND 1000
#define TOTAL_COLORS 3 //total 3 colors, RGB

Adafruit_TCS34725 tcs = Adafruit_TCS34725(TCS34725_INTEGRATIONTIME_600MS, TCS34725_GAIN_1X);

unsigned long startTime_avgVal = 0; //start time
unsigned long startTime_totalVal = 0; //start time
float totalVal = 0;
int counter = 0;

void setup(void) {
  Serial.begin(9600);
  // pinMode(ONBOARD_LED,OUTPUT);
  pinMode(RED_LED,OUTPUT);
  pinMode(YELLOW_LED,OUTPUT);
  pinMode(GREEN_LED,OUTPUT);

  if (tcs.begin()) {
    Serial.println("Found sensor");
  } else {
    Serial.println("No TCS34725 found ... check your connections");
    while (1);
  }

}

void loop(void) {
  float r, g, b, avgVal,totalAvgValue;
  unsigned long currentTime = millis();

  if(currentTime - startTime_avgVal > ONE_SECOND){
    timePassed(); Serial.println("One second completed");
    tcs.getRGB(&r, &g, &b);
    
    //reading sensor values
    timePassed(); Serial.print("R: "); Serial.println(r);
    timePassed(); Serial.print("G: "); Serial.println(g);
    timePassed(); Serial.print("B: "); Serial.println(b);
    avgVal = (r+g+b)/TOTAL_COLORS; //average values
    timePassed(); Serial.print("Average value: "); Serial.println(avgVal);
    totalVal = totalVal + avgVal;
    counter++; //number of readings
    startTime_avgVal = currentTime;
  }
  //checking for two seconds
  if(currentTime - startTime_totalVal > TWO_SECONDS){
    timePassed(); Serial.println("Two second completed");
    timePassed(); Serial.print("counter: "); Serial.println(counter);
    timePassed(); Serial.print("total value: "); Serial.println(totalVal);
    timePassed(); Serial.print("totalAvgValues"); Serial.println(totalVal/counter);
    timePassed(); Serial.println("Two seconds, yayyy!! 1");
    totalAvgValue = (totalVal/counter);
    if(totalAvgValue>MAX_THRESHOLD){
      //RED
      timePassed(); Serial.println("RED LED ON");
      digitalWrite(RED_LED,HIGH);
      digitalWrite(GREEN_LED,LOW);
      digitalWrite(YELLOW_LED,LOW);
    }else if(MIN_THRESHOLD < totalAvgValue ){
      //YELLOW
      timePassed(); Serial.println("YELLOW LED ON");
      digitalWrite(YELLOW_LED,HIGH);
      digitalWrite(GREEN_LED,LOW);
      digitalWrite(RED_LED,LOW);      
    }else if(totalAvgValue<MIN_THRESHOLD){
      //GREEN
      timePassed(); Serial.println("GREEN LED ON");
      digitalWrite(GREEN_LED,HIGH);
      digitalWrite(RED_LED,LOW);
      digitalWrite(YELLOW_LED,LOW);
    }
  startTime_totalVal = currentTime;
  counter=0;
  totalVal=0;
  }
  timePassed(); Serial.println(" ");
}

/*
to print time passed after reboot in ms
*/
void timePassed(){
  // return millis();
  Serial.print("Time after reboot in ms: "); Serial.print(millis()); Serial.print(" >> ");
}