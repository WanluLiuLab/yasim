
export PROJECT_NAME := proc_profiler

include .maint/makefiles/common_environment.mk
include .maint/makefiles/general_pub.mk
include .maint/makefiles/python_test.mk
include .maint/makefiles/python_pub.mk
include .maint/makefiles/doc_adapter.mk

clean:
	$(RM) .pytest_cache

distclean: clean
	$(RM) activate.sh renv.lock .Rprofile renv venv

