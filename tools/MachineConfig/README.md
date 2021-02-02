# MachineConfig

Скрипт для настройки машины (биоинформатика, офис).

### Дополнительная настройка

#### ffmpeg

Скрипт для упрощения нарезки видео:

```bash
#!/bin/bash

IN_FILE="The.100.S03E04.720p.WEB-DL.Rus.Eng.HDCLUB.mkv"; OUT_FILE="/_home/Temp/roan.mkv"; TIME_START="00:27:25"; TIME_END="00:31:31"; V_TRACK=0; A_TRACK=2; ffmpeg -i "$IN_FILE" -qscale 0 -filter_complex "[0:v:0]subtitles="$IN_FILE":si=0[v]" -map "[v]" -map 0:a:$A_TRACK -ss $TIME_START -to $TIME_END "$OUT_FILE"
```
-f lavfi -i color=c=black:s=1920x1080 -filter_complex "[0:v]scale=w=iw:h=ih[scaled];[1:v][scaled]overlay=x=0.10*main_w:y=0.10*main_h:eof_action=endall[out];[0:a]anull[aud]" -map "[out]" -map "[aud]" -strict -2
