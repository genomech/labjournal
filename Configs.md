# Конфигурация машины

Ubuntu, разумеется.

## privoxy + tor

Связка прокси-сервера с выходом в луковую сеть.

Установка:

```bash
sudo apt install tor tor-geoipdb privoxy
```

Если не устанавливается privoxy, нужно следовать указаниям на [официальном сайте](https://www.privoxy.org/).

Если не устанавливается tor:

```bash
ver=$(lsb_release -c -s); echo "deb https://deb.torproject.org/torproject.org "$ver" main" | sudo tee -a /etc/apt/sources.list

curl https://deb.torproject.org/torproject.org/A3C4F0F979CAA22CDBA8F512EE8CBC9E886DDD89.asc | gpg --import

gpg --export A3C4F0F979CAA22CDBA8F512EE8CBC9E886DDD89 | sudo apt-key add -

sudo apt update

sudo apt install tor tor-geoipdb
```

Настройка:

```bash
sudo mv /etc/privoxy/config /etc/privoxy/config.backup

echo -e "forward-socks5 / 127.0.0.1:9050 .\nconfdir /etc/privoxy\nlogdir /var/log/privoxy\nactionsfile default.action\nactionsfile user.action\nfilterfile default.filter\nlogfile logfile\ndebug 4096\ndebug 8192\nuser-manual /usr/share/doc/privoxy/user-manual\nlisten-address 127.0.0.1:8118\ntoggle 1\nenable-remote-toggle 0\nenable-edit-actions 0\nenable-remote-http-toggle 0\nbuffer-limit 4096" | sudo tee /etc/privoxy/config

sudo service privoxy restart

sudo service tor restart
```

Далее можно использовать как прокси-сервер 127.0.0.1:9050.
Работает как с торрент-клиентами, так и с Telegram.
