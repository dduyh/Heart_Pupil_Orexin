// triggered by Ethovision output, then start 40hz pulse (10+15) or 20hz(15+35)
int triggerPin =2;
byte LaserPIN = 3;
void setup() {
  pinMode(triggerPin, INPUT);
  pinMode(LaserPIN, OUTPUT);
}
void loop(){
 if (digitalRead(triggerPin)==HIGH){
   unsigned long OFFtimings = 45; // milliseconds 35 is 20hz
   unsigned long ONTimings = 5; // milliseconds
   digitalWrite(LaserPIN, HIGH);
   delay(ONTimings);
   digitalWrite(LaserPIN, LOW);
   delay(OFFtimings);
  }  //End of if()
} //End of loop()
