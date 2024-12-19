int triggerPin = 2;
byte LaserPIN = 3;

// Define ON and OFF times
unsigned long OFFtimings = 45; // milliseconds, 45ms for 20Hz
unsigned long ONtimings = 5;  // milliseconds

void setup() {
  pinMode(triggerPin, INPUT);
  pinMode(LaserPIN, OUTPUT);
}

void loop() {
  // Check if TTL signal is HIGH
  if (digitalRead(triggerPin) == HIGH) {
    // Output 20Hz pulses for 10 seconds
    unsigned long startTime = millis();  // Record the start time
    while (millis() - startTime < 10000) {  // 10 seconds
      digitalWrite(LaserPIN, HIGH);
      delay(ONtimings);
      digitalWrite(LaserPIN, LOW);
      delay(OFFtimings);

      // Check if TTL signal has changed to LOW, exit loop if it has
      if (digitalRead(triggerPin) == LOW) {
        break;
      }
    }

    // Turn off laser for 20 seconds
    digitalWrite(LaserPIN, LOW);  // Ensure laser is off
    delay(20000);
  } else {
    // Ensure laser is off when TTL is LOW
    digitalWrite(LaserPIN, LOW);
  }
}
