from spack import *
import os

class Quicksched(BundlePackage):
    """dummy placeholder for oomph dependencies"""
    homepage    = "https://127.0.0.1/readme.html"
    git = "https://github.com/biddisco/QuickSched"
    generator   = "Ninja"
    maintainers = ["biddisco"]

    version("master")

    # ------------------------------------------------------------------------
    # variants
    # ------------------------------------------------------------------------
    #variant("mpi", default=False, description="Enable mpi for pika dependency")

    # ------------------------------------------------------------------------
    # build time dependencies
    # ------------------------------------------------------------------------
    depends_on("ninja", type="build")
    depends_on("cmake@3.23:", type="build")

    # ------------------------------------------------------------------------
    # generic c++ libs
    # ------------------------------------------------------------------------
    depends_on("boost +atomic+chrono+container+context+coroutine+date_time+filesystem+program_options+regex+serialization+system+test+thread+json cxxstd=17")
    depends_on("stdexec@main")
    depends_on("fmt")

    # ------------------------------------------------------------------------
    # IO
    # ------------------------------------------------------------------------
    depends_on("mpi")
    depends_on("hdf5 +mpi")

    # ------------------------------------------------------------------------
    # maths
    # ------------------------------------------------------------------------
    depends_on("openblas") 

    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    def cmake_args(self):
        spec = self.spec

        return []
