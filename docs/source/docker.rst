Docker Image
============

The riptide `Dockerfile`_ is located in the ``docker`` subdirectory. To build the image, clone the repository and in its base directory type:

.. _`Dockerfile`: https://github.com/v-morello/riptide/blob/master/docker/Dockerfile

.. code-block:: console

    make docker

Which builds an image named ``riptide-ffa``. Both python and ipython are installed within the docker image. To start a container:

.. code-block:: console

    docker run -it --rm riptide-ffa

Feel free to adapt the Dockerfile to your needs. Refer to your favourite docker cheat sheet for further information.