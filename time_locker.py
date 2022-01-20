import time
from ctypes import *
from ctypes.wintypes import *

my_minute_1 = '29'# 设置时间，可以设定在多个时间点锁屏，下面的判断条件改一下就行
my_minute_2 = '59'

def main():
  shell32 = windll.LoadLibrary("shell32.dll")
  while True:
    t = time.localtime() # 当前时间的纪元值
    minute = time.strftime("%M", t) # 将纪元值转化为包含时、分的字符串
    if minute == my_minute_1 or minute == my_minute_2:
      shell32.ShellExecuteW(None,'open', 'rundll32.exe','USER32,LockWorkStation','',5)#调用系统锁屏
      #如果不想强制锁定，只锁定一次的话，把下面这句加上就行
      #time.sleep(60)
    time.sleep(1)#暂停一秒，节省资源
 
if __name__ == "__main__":
  print("程序将在每小时 {} 分和 {} 锁定系统，起来走动一下，喝点水哦".format(my_minute_1, my_minute_2))
  main()