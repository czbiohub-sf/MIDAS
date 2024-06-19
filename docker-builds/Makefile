NAME=midas
VERSION:=$(shell cat ../midas/__init__.py  | grep "^version" | awk -F '\"' '{print "v"v$2}')

REGISTRY=zhaoc1/$(NAME)


all: build push run

build:
    docker image build -t $(NAME):$(VERSION) -t $(NAME):latest -f Dockerfile .
    docker image tag $(NAME):latest $(REGISTRY):$(VERSION)
    docker image tag $(NAME):latest $(REGISTRY):latest

push:
    docker image push $(REGISTRY):$(VERSION)
    docker image push $(REGISTRY):latest

run:
    docker container run --rm $(NAME):latest
    docker run --volume "/home/ubuntu/.aws":"/root/.aws":ro --rm -it $(NAME):latest
