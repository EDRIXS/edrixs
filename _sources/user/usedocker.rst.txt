.. _edrixsanddocker:

*********************
edrixs and docker
*********************

Run edrixs in a docker container
--------------------------------

To make life easier, we have built a docker image based on Ubuntu Linux (22.04) for edrixs, so you don't need to struggle with the installation anymore.
The docker image can be used on any OS as long as the `docker <https://www.docker.com/>`_ application is available.
Follow these steps to use the docker image:

* Install the `docker <https://www.docker.com/>`_ application on your system.

* Once Docker is running, create a directory to store data and create a file called ``docker-compose.yml`` with contents ::

    version:  '3'
    services:
      edrixs-jupyter:
          image: edrixs/edrixs
          volumes:
            - ./:/home/rixs
          working_dir: /home/rixs
          ports:
            - 8888:8888
          command: "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root"
      edrixs-ipython:
          image: edrixs/edrixs
          volumes:
            - ./:/home/rixs
          working_dir: /home/rixs

  and execute ::

    docker compose up

  This will return a url, which you can open to connect to the jupyter session. 

* If you would like to access a terminal rather than jupyter run ::

    docker compose run --rm edrixs-ipython


Sharing your code
-----------------

Using Docker is a nice way to straightforwardly share your code with others. The standard way to specify which docker image is needed to run your code is to include a file named ``Dockerfile`` with the following contents ::

    FROM edrixs/edrixs

You might like to checkout the `jupyter-repo2docker
<https://repo2docker.readthedocs.io/en/latest/>`_ project, which helps automate the process of building and connecting to docker images. The `mybinder <https://mybinder.org/>`_ project might also be helpful as this will open a github respository of notebooks in an executable environment, making your code immediately reproducible by anyone, anywhere.
