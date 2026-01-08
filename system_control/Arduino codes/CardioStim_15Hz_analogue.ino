byte LaserPIN = 3;
bool isOn = false;
char terminator = '/';  // Set terminator for each input command
String mode;            // mode string saves the command mode

// Define ON and OFF times
unsigned long OFFtimings = 57; // milliseconds 57 is 15hz
unsigned long ONTimings = 10; // milliseconds
unsigned long lastTime = 0;
bool ledState = false;

void setup() {

  Serial.begin(9600);  // Begin serial communication (baud rate 9600) with MATLAB

  pinMode(LaserPIN, OUTPUT);

  analogWrite(LaserPIN, 0);  // Initiate the output for LED
}

void loop() {

  if (Serial.available() > 0) {  // Wait for user input

    mode = Serial.readStringUntil(terminator);  // Read the input command from MATLAB to set mode

    if (mode == "led_on")  // The mode is turning on the LEDs
    {
      isOn = true;
      lastTime = millis();
      ledState = false;
      analogWrite(LaserPIN, 0);
    }
    else if (mode == "led_off") {    // The mode is turning off the LEDs
      isOn = false;
      analogWrite(LaserPIN, 0);
    }
  }


  if (isOn) {
    unsigned long now = millis();
    if (!ledState && now-lastTime>=OFFtimings) {
      analogWrite(LaserPIN, 168);
      ledState = true;
      lastTime = now;
    }
    else if (ledState && now-lastTime>=ONTimings) {
      analogWrite(LaserPIN, 0);
      ledState = false;
      lastTime = now;   
    }
  }
} //End of loop()
