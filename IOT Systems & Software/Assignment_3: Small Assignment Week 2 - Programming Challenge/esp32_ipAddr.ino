#include <LiquidCrystal_I2C.h>
#include <WiFi.h>

// set the LCD number of columns and rows
int lcdColumns = 16;
int lcdRows = 2;
const char* ssid = "ResNet Mobile Access";
#define ONBOARD_LED 2

// set LCD address, number of columns and rows
LiquidCrystal_I2C lcd(0x27, lcdColumns, lcdRows);  

void setup(){
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid);

  pinMode(ONBOARD_LED,OUTPUT); //settings for led

  Serial.print("Connecting to WiFi..");
  while(WiFi.status() != WL_CONNECTED){
    digitalWrite(ONBOARD_LED,LOW);
    delay(10);
    Serial.print('.');
  }
  digitalWrite(ONBOARD_LED,HIGH);

  Serial.println("");
  Serial.println("Connected to WiFi network!");
  Serial.print("IP address: ");
  Serial.println(WiFi.localIP());
  // initialize LCD
  lcd.init();
  // turn on LCD backlight                      
  lcd.backlight();
  lcd.setCursor(0, 0);
  lcd.print("IP Address");
}

void loop(){
  lcd.setCursor(0, 1);
  lcd.print(WiFi.localIP());
}