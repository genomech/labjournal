# MachineConfig

Скрипт для настройки машины (биоинформатика, офис).

### Дополнительная настройка

#### ffmpeg

Скрипт для упрощения нарезки видео:

```bash
#!/bin/bash

IN_FILE="/media/avicenna/Dopamine_Python/GoT/GOT/Season 08/GOT.[S08E02].2xRu.En.[qqss44].mkv"
OUT_FILE="/dev/my/MyDocs/temp/jenny.mkv"
TIME_START="00:51:05"
TIME_END="00:53:01"
V_TRACK=0
A_TRACK=2

ffmpeg -i "$IN_FILE" -qscale 0 -map 0:v:$V_TRACK -map 0:a:$A_TRACK -ss $TIME_START -to $TIME_END "$OUT_FILE"
```
