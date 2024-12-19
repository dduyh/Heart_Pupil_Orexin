/********************************************************************************************
*        Date:    November 30th, 2023                                                       *
*       Author:  Yihui Du                                                                   *
*  Description:  valve_driver                                                               *
*                                                                                           *
********************************************************************************************/

int valve = 2;  // Set the digital pin 2 to control valve driver

char terminator = '/';  // Set terminator for each input command
String mode;            // mode string saves the command mode

// Setup code are put here to run once:
void setup() {

  Serial.begin(9600);  // Begin serial communication (baud rate 9600) with MATLAB

  pinMode(valve, OUTPUT);  // sets the digital pin 2 as output

  digitalWrite(valve, LOW);  // Initiate the output for valve
}

// Main code are put here to run repeatedly:
void loop() {
  while (Serial.available() == 0) {  // Wait for user input
  }
  mode = Serial.readStringUntil(terminator);  // Read the input command from MATLAB to set mode

  if (mode == "first_on")  // The mode is turning on the first valve drivers
  {
    digitalWrite(valve, HIGH);  // Set the digital left pin to HIGH
  }
  if (mode == "first_off")  // The mode is turning off the first valve drivers
  {
    digitalWrite(valve, LOW);   // Set the digital left pin to LOW
  }
}
