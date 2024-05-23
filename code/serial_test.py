import serial
import time

ser = serial.Serial(port='/dev/cu.usbserial-56470095751')
# colours = ["violet", "royal_blue", "blue", "cyan", "green", "lime", "amber", "red_orange", "red", "deep_red", "far_red", "white", "mint"]
# read=[]
# while True:
# 	r = ser.read_until()
# 	clean_r = int(str(r)[2:-3])
# 	read.append(clean_r)
# 	if len(read) > 13:
# 		break
# 	print(clean_r)
# 	print(colours[clean_r])

# for i in range(13):
#     ser.write(i)
#     time.sleep(2)
    
while True:
    input_value = input('Enter color: ')
    ser.write(input_value.encode())
    # print("write: ", input_value.encode())
    # print('readed:', ser.read(3))