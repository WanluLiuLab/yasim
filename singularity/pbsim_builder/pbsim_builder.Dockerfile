FROM alpine:latest
RUN apk update
RUN apk add gcc g++ binutils libstdc++-dev libstdc++
