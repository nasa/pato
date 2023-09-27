TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -DOPENMP -fopenmp -lpthread
QMAKE_CXXFLAGS += -DGL_GLEXT_PROTOTYPES  -DOPENMP -fopenmp --std=c++0x
QMAKE_LFLAGS +=  -fopenmp

SOURCES += \
    ../testframework/*.cpp\
    ../testsuites/*.cpp \
    ../*.cpp

HEADERS += \
    ../testframework/*.h \
    ../*.h


