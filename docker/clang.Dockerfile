FROM alpine:edge

MAINTAINER Rajeeva Musunuri <rmusunuri@nygenome.org>

RUN apk update && apk upgrade && apk add --no-cache ca-certificates && rm -rf /var/cache/apk/* && \
    apk add --no-cache bash make cmake ninja git linux-headers libc-dev gcc g++ clang \
                       zlib-dev bzip2-dev xz-dev curl-dev openssl-dev

ENV CC  clang
ENV CXX clang++

RUN mkdir -p /usr/src && git clone https://github.com/omicsnut/v2_lancet.git /usr/src/lancet && \
    mkdir -p /usr/src/lancet/build && cd /usr/src/lancet/build && \
    cmake -GNinja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release .. && \
    ninja -v && cd / && cp /usr/src/lancet/build/lancet /usr/bin/ && find /usr/src \! -name "*mimalloc-build*" -delete

ENTRYPOINT ["lancet"]
