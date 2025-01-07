import time
import pyautogui as pg
import pyperclip as pc
#from apscheduler.schedulers.blocking import BlockingScheduler

# 操作间隔为1秒
pg.PAUSE = 1

name = '小玲犟'
msg = '给老子，我看看能发送不'
send_time = '15:58:00'


def main():
    # 打开微信
    pg.hotkey('ctrl', 'alt', 'w')
    pg.hotkey('ctrl', 'f')

    # 找到女朋友
    pc.copy(name)
    pg.hotkey('ctrl', 'v')
    pg.press('enter')

    # 发送消息
    pc.copy(msg)
    pg.hotkey('ctrl', 'v')
    pg.press('enter')

    # 隐藏微信
    pg.hotkey('ctrl', 'alt', 'w')


if __name__ == '__main__':
    # 循环获取当前时间，如果为指定时间我们将调用main函数为女朋友发送晚安消息
    while True:
        hour = time.localtime()
        now_time = time.strftime("%H:%M:%S", hour)
        # print(f'\r{now_time}')
        # 如果时间为半夜12点，则给女朋友发消息
        if now_time == send_time:
            print(f'时间已到：{send_time}，开始为女朋友发送消息！')
            main()
            print("任务结束！")
            break