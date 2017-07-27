avogadro_plugin(Ribbon
  "Ribbon rendering scheme"
  ScenePlugin
  ribbon.h
  Ribbon
  ribbon.cpp
  "")

target_link_libraries(Ribbon LINK_PRIVATE AvogadroRendering)
