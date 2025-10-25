#include <WiFi.h>
#include <WiFiUdp.h>

WiFiUDP udp;

char packetBuffer[255];
unsigned int localPort = 9999;                
const char *serverip = "192.168.137.101"; // Raspberry Pi IP
unsigned int serverport = 8888;// Pi's UDP port

const char *ssid = "deepanshu_laptop";        
const char *password = "anshu.com";           

const int LED_PIN = 2;// Onboard LED
const int LIGHT_SENSOR_PIN = 34;              

// Timing variables
unsigned long previousBlink = 0;
unsigned long previousLight = 0;
unsigned long previousSend = 0;

const long blinkInterval = 500;// LED blink every 0.5s
const long lightInterval = 1000;// LDR sample every 1s
const long sendInterval = 2000;// Send avg every 2s

#define NUM_SAMPLES 5
int samples[NUM_SAMPLES] = {0};
int sampleIndex = 0;
int sampleCount = 0;

bool isRunning = false;   // false = idle, true = running

//setup
void setup() {
  Serial.begin(115200);
  pinMode(LED_PIN, OUTPUT);
  digitalWrite(LED_PIN, LOW);

  Serial.println("\nConnecting to WiFi...");
  WiFi.begin(ssid, password);
  while (WiFi.status() != WL_CONNECTED) {
    delay(500);
    Serial.print(".");
  }

  Serial.println("\nWiFi connected!");
  Serial.print("ESP32 IP: ");
  Serial.println(WiFi.localIP());

  udp.begin(localPort);
  Serial.printf("UDP listening on port %d\n", localPort);
  Serial.println("Waiting for command from Raspberry Pi...");
}

//reset state
void resetState() {
  isRunning = false;
  digitalWrite(LED_PIN, LOW);// Turn off LED

  previousBlink = millis();
  previousLight = millis();
  previousSend = millis();

  sampleIndex = 0;
  sampleCount = 0;
  for (int i = 0; i < NUM_SAMPLES; i++) samples[i] = 0;

  Serial.println("[ESP] Reset to initial state, waiting for new command...");
}

//main loop
void loop() {
  unsigned long currentMillis = millis();

  //checking for incoming UDP
  int packetSize = udp.parsePacket();
  if (packetSize) {
    int len = udp.read(packetBuffer, 255);
    if (len > 0) packetBuffer[len] = '\0';
    Serial.print("[Received from Pi] ");
    Serial.println(packetBuffer);

    if (strstr(packetBuffer, "Start Communication") != NULL) {
      if (!isRunning) {
        isRunning = true;
        previousBlink = millis();
        previousLight = millis();
        previousSend = millis();
        Serial.println("[ESP] Starting sensor reading and LED blink...");
      }
    } else if (strstr(packetBuffer, "Reset/Error Cleared") != NULL) {
      if (isRunning) {
        resetState();// Stop tasks and reset
      }
    } else {
      Serial.println("[ESP] Unknown command, ignoring.");
    }
  }
  if (!isRunning) return;

  //blink LED 0.5 sec
  if (currentMillis - previousBlink >= blinkInterval) {
    previousBlink = currentMillis;
    digitalWrite(LED_PIN, !digitalRead(LED_PIN));
  }

  //Read LDR every 1 sec
  if (currentMillis - previousLight >= lightInterval) {
    previousLight = currentMillis;

    int sensorValue = analogRead(LIGHT_SENSOR_PIN);
    samples[sampleIndex] = sensorValue;
    sampleIndex = (sampleIndex + 1) % NUM_SAMPLES;
    if (sampleCount < NUM_SAMPLES) sampleCount++;

    Serial.printf("[Light Sensor] Value: %d\n", sensorValue);
  }

  //Send average every 2 sec
  if (sampleCount == NUM_SAMPLES && (currentMillis - previousSend >= sendInterval)) {
    previousSend = currentMillis;

    long sum = 0;
    for (int i = 0; i < NUM_SAMPLES; i++) sum += samples[i];
    float avg = sum / (float)NUM_SAMPLES;

    Serial.printf("[Average] Last %ds Avg = %d\n", NUM_SAMPLES, (int)avg);

    //Send average via UDP
    char msg[100];
    snprintf(msg, sizeof(msg), "Average Light=%d", (int)avg);
    udp.beginPacket(serverip, serverport);
    udp.print(msg);
    udp.endPacket();
  }
}
