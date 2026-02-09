************************
Developing within docker
************************
To make changes to edrixs, it can be useful to work within a docker environment. Here are some instructions for doing so.

* Fork edrixs to your GitHub account and
  and execute ::

    git clone https://github.com/<YOUR_USERNAME>/edrixs.git
    cd edrixs
    git checkout -b <YOUR_NEW_BRANCH>

* Build a docker image from within the edrixs directory targeting the developer stage and tagging the image as ``edrixs_developer`` ::

    docker build --target developer -f docker/Dockerfile -t edrixs_developer .

* Go one level up in the directory structure and create a file called ``docker-compose.yml`` with contents ::

    version:  '3'
    services:
    edrixs-jupyter:
        image: edrixs_developer
        volumes:
            - ./:/home/rixs
        working_dir: /home/rixs
        ports:
            - 8888:8888
        command: "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root"
    edrixs-ipython:
        image: edrixs_developer
        volumes:
            - ./:/home/rixs
        working_dir: /home/rixs

  Executing ::

    docker compose up

  will generate a jupyter session. Executing ::

    docker compose run --rm edrixs-ipython

  will generate an interactive terminal session.

* From within the container, reinstalling edrixs, running the tests, and building the documentation can be done via the following commands::

    pip install -v .

    ./scripts/run-tests.sh

    make -C docs html
