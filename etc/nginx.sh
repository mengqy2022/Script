#!/bin/bash

# 定义Nginx服务的相关路径
NGINX_CONF="/etc/nginx/nginx.conf"
ACCESS_LOG="/var/log/nginx/access.log"
ERROR_LOG="/var/log/nginx/error.log"

# 执行Nginx服务操作
nginx_operation() {
    local action=$1
    if sudo service nginx $action; then
        echo "Nginx服务已$action"
    else
        echo "执行Nginx服务$action时出错"
    fi
}

# 检查Nginx配置文件语法
check_nginx_config() {
    local result
    result=$(sudo nginx -t 2>&1)
    if echo "$result" | grep -q "syntax is ok"; then
        echo "Nginx配置文件语法正确"
    else
        echo "Nginx配置文件存在语法错误，请检查"
        echo "$result"
    fi
}

# 备份Nginx配置文件
backup_nginx_config() {
    local backup_file="${NGINX_CONF}.$(date +%Y%m%d%H%M%S).bak"
    if sudo cp $NGINX_CONF $backup_file; then
        echo "Nginx配置文件已备份至 $backup_file"
    else
        echo "备份Nginx配置文件时出错"
    fi
}

# 重新加载Nginx配置
reload_nginx_config() {
    if sudo nginx -s reload; then
        echo "Nginx配置已重新加载"
    else
        echo "重新加载Nginx配置时出错"
    fi
}

# 查看Nginx日志
view_log() {
    local log_file=$1
    if [ -f "$log_file" ]; then
        sudo tail -f $log_file
    else
        echo "日志文件 $log_file 不存在"
    fi
}

# 清理Nginx访问日志
clean_access_log() {
    if sudo truncate -s 0 $ACCESS_LOG; then
        echo "Nginx访问日志已清理"
    else
        echo "清理Nginx访问日志时出错"
    fi
}

# 检查Nginx是否存在安全漏洞
check_nginx_security() {
    local nginx_ip
    read -p "请输入Nginx服务器的IP地址: " nginx_ip
    if ! command -v nmap &> /dev/null; then
        echo "nmap未安装，请先安装nmap"
    elif ! nmap -sV --script http-nginx-version $nginx_ip; then
        echo "nmap执行失败，请检查nmap是否已安装并正确配置"
    fi
}

# 设置Nginx的访问控制
set_nginx_access_control() {
    local deny_ip
    read -p "请输入要拒绝的IP地址或IP地址段（格式：IP 或 IP/掩码）: " deny_ip
    if ! grep -q "deny $deny_ip;" $NGINX_CONF; then
        if sudo sed -i "/server {/a\    deny $deny_ip; allow all;" $NGINX_CONF; then
            echo "已设置拒绝 $deny_ip 访问，允许其他所有IP访问"
        else
            echo "设置Nginx访问控制时出错"
        fi
    else
        echo "该IP地址或IP地址段已存在于Nginx配置中，无需重复添加"
    fi
}

# 主菜单
display_menu() {
    echo "请选择要执行的操作："
    echo "1. 启动Nginx服务"
    echo "2. 停止Nginx服务"
    echo "3. 重启Nginx服务"
    echo "4. 检查Nginx配置文件语法"
    echo "5. 备份Nginx配置文件"
    echo "6. 加载新的Nginx配置"
    echo "7. 查看Nginx访问日志"
    echo "8. 清理Nginx访问日志"
    echo "9. 查看Nginx错误日志"
    echo "10. 添加后端服务器到负载均衡池"
    echo "11. 从负载均衡池中移除后端服务器"
    echo "12. 重新加载负载均衡配置"
    echo "13. 检查Nginx是否存在安全漏洞"
    echo "14. 设置Nginx的访问控制"
    echo "15. 退出"
    read -p "请输入选项编号: " choice
}

# 主程序逻辑
while true; do
    display_menu
    case $choice in
        1) nginx_operation start ;;
        2) nginx_operation stop ;;
        3) nginx_operation restart ;;
        4) check_nginx_config ;;
        5) backup_nginx_config ;;
        6) reload_nginx_config ;;
        7) view_log $ACCESS_LOG ;;
        8) clean_access_log ;;
        9) view_log $ERROR_LOG ;;
        10) add_backend_server ;;
        11) remove_backend_server ;;
        12) reload_nginx_config ;;
        13) check_nginx_security ;;
        14) set_nginx_access_control ;;
        15) echo "感谢使用，再见！"; break ;;
        *) echo "无效的选项，请重新输入。" ;;
    esac
done
