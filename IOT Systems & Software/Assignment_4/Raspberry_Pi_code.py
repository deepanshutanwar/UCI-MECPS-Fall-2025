import socket
import threading
from gpiozero import LED, Button
import time

LOCAL_PORT = 8888
ESP_IP = "192.168.137.219"  #ESP32 IP
ESP_PORT = 9999
BUFFER_SIZE = 1024
TIMEOUT_SEC = 10
BLINK_INTERVAL = 0.5

BUTTON_PIN = 17
LED_PIN = 18
led = LED(LED_PIN)
button = Button(BUTTON_PIN, pull_up=True)

#RGB LEDs
RED_PIN = 23
GREEN_PIN = 24
BLUE_PIN = 16
red_led = LED(RED_PIN)
green_led = LED(GREEN_PIN)
blue_led = LED(BLUE_PIN)

button_pressed = False
error_mode = False
last_event_time = None
lock = threading.Lock()

#setup UDP
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

def init_udp():
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
    sock.bind(('', LOCAL_PORT))
    print(f"UDP client started on port {LOCAL_PORT}, target ESP32: {ESP_IP}:{ESP_PORT}")

#LED control
def led_control():
    while True:
        with lock:
            if error_mode:
                # Blink main LED, RGB LEDs off
                led.toggle()
                red_led.off()
                green_led.off()
                blue_led.off()
            elif button_pressed:
                led.on()      # Steady ON after first button press
            else:
                led.off()     # Default OFF
        time.sleep(BLINK_INTERVAL)

#watchdog
def watchdog():
    global error_mode
    while True:
        time.sleep(0.1)
        with lock:
            if button_pressed and not error_mode and last_event_time:
                elapsed = time.time() - last_event_time
                if elapsed >= TIMEOUT_SEC:
                    print(f"[ERROR] No response from ESP32 for {TIMEOUT_SEC} seconds ? LED blinking")
                    error_mode = True  # Start blinking

#update RGB LEDs
def update_rgb_leds(light_value: int):
    if error_mode:
        return

    red_led.off()
    green_led.off()
    blue_led.off()

    if light_value < 1000:
        red_led.on()
        print("RGB LEDs ? LOW (Red ON)")
    elif 1000 <= light_value <= 2000:
        red_led.on()
        green_led.on()
        print("RGB LEDs ? MEDIUM (Red + Green ON)")
    else:  # >2000
        red_led.on()
        green_led.on()
        blue_led.on()
        print("RGB LEDs ? HIGH (All ON)")



#UDP Receiver
def receiver():
    global last_event_time, error_mode
    while True:
        data, addr = sock.recvfrom(BUFFER_SIZE)
        msg = data.decode().strip()
        print(f"\n[From {addr[0]}:{addr[1]}] {msg}")
        with lock:
            last_event_time = time.time()
            # Reset error mode
            if error_mode:
                print("[INFO] Data received from ESP32 ? stopping main LED blink")
                error_mode = False

            # Parsing average value
            if "Average Light=" in msg:
                try:
                    light_val = int(msg.split("=")[1])
                    print(f"Light sensor value received: {light_val}")
                    update_rgb_leds(light_val)
                except ValueError:
                    print("[ERROR] Invalid light sensor value received")

#Button handler
def on_button_press():
    global button_pressed, last_event_time, error_mode
    with lock:
        if error_mode:
            msg = "Button Pressed - Reset/Error Cleared"
            print("\n[Button Pressed] Clearing error ? sending reset message to ESP and turning LED OFF")
            error_mode = False  # Stop blinking
            button_pressed = False  # LED OFF after reset
            led.off()
        else:
            msg = "Button Pressed - Start Communication"
            print("\n[Button Pressed] Starting communication ? sending message to ESP")
            button_pressed = True
            led.on()  # Steady ON after first press

        sock.sendto(msg.encode(), (ESP_IP, ESP_PORT))
        print(f"Sent to {ESP_IP}:{ESP_PORT} ? {msg}")
        last_event_time = time.time()  # Reset timer

#Main
if __name__ == "__main__":
    try:
        init_udp()

        threading.Thread(target=receiver, daemon=True).start()
        threading.Thread(target=watchdog, daemon=True).start()
        threading.Thread(target=led_control, daemon=True).start()

        button.when_pressed = on_button_press

        #optional manual
        while True:
            msg = input("You: ")
            if msg.strip():
                sock.sendto(msg.encode(), (ESP_IP, ESP_PORT))

    except KeyboardInterrupt:
        print("\nExiting gracefully...")
    finally:
        sock.close()
        led.off()

