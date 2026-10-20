#!/bin/bash
### Post a message and some files to Slack.
### Usage:  slack_upload.sh "message" file1 [file2 ...]
### Files that don't exist or fail to upload are listed in the message.

SLACK_CHANNEL=C0CB96QL6TA                  # channel ID (not secret)
SLACK_TOKEN_FILE=~/.slack_token            # bot token (xoxb-...), readable only by you

MSG=$1; shift

### The auth header goes in a private temp file, so the token never appears on a
### command line where other users could see it with ps
SLACK_AUTH=$(mktemp)
trap 'rm -f "$SLACK_AUTH"' EXIT
printf 'Authorization: Bearer %s\n' "$(< $SLACK_TOKEN_FILE)" > "$SLACK_AUTH"

### Upload each file (Slack's two-step upload: get an upload URL, then send the file)
FILE_IDS=""; MISSING=""; FAILED=""
for f in "$@"; do
    [ -f "$f" ] || { MISSING="$MISSING $(basename "$f")"; continue; }
    read -r URL FID < <(curl -s -H @"$SLACK_AUTH" https://slack.com/api/files.getUploadURLExternal \
                             -d filename="$(basename "$f")" -d length=$(stat -c %s "$f") |
                        python3 -c 'import json,sys; d=json.load(sys.stdin); print(d.get("upload_url","-"), d.get("file_id",""))' 2>/dev/null)
    if [[ $URL == http* ]] && curl -sf -o /dev/null -F filename=@"$f" "$URL"; then
        FILE_IDS="$FILE_IDS{\"id\":\"$FID\"},"
    else
        FAILED="$FAILED $(basename "$f")"
    fi
done

### Share the uploaded files in the channel with one message (just a message if none uploaded)
[ -n "$MISSING" ] && MSG="$MSG. Not produced:$MISSING"
[ -n "$FAILED" ]  && MSG="$MSG. Upload failed:$FAILED"
if [ -n "$FILE_IDS" ]; then
    SLACK_METHOD=files.completeUploadExternal
    SLACK_BODY="{\"files\":[${FILE_IDS%,}],\"channel_id\":\"$SLACK_CHANNEL\",\"initial_comment\":\"$MSG\"}"
else
    SLACK_METHOD=chat.postMessage
    SLACK_BODY="{\"channel\":\"$SLACK_CHANNEL\",\"text\":\"$MSG\"}"
fi
curl -s -H @"$SLACK_AUTH" -H 'Content-Type: application/json; charset=utf-8' \
     https://slack.com/api/$SLACK_METHOD -d "$SLACK_BODY" |
    python3 -c 'import json,sys; d=json.load(sys.stdin); print("Slack upload:", "ok" if d.get("ok") else d.get("error"))' 2>/dev/null ||
    echo "Slack upload: no response from Slack"
