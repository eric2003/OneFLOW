Tecplot
==================================

#. `Tecplot hoempage <https://www.tecplot.com/>`_
#. `Source code for tecplot.data.dataset <https://www.tecplot.com/docs/pytecplotdocs/v0.7/_modules/tecplot/data/dataset.html>`_
#. `Tecplot 360 Data Format Guide <http://www.wagnerrp.com/files/plugin-dataformat.pdf>`_
#. `Higher-Order Elements in Tecplot 360 Beta <https://www.youtube.com/watch?v=koK6rvJAbto/>`_
#. `foamToTecplot360 <https://searchcode.com/codesearch/view/65520755/>`_
#. `mat2tecplot.m <https://github.com/wme7/aero-matlab/blob/master/Tecplot/mat2tecplot.m>`_
#. `TP - TecPlot ascii file format <https://paulbourke.net/dataformats/tp/>`_
#. `Using MATLAB® & TecIO to Read/Write Tecplot Data File Formats <https://www.youtube.com/watch?v=PofNOjBK7Z8/>`_
#. `Back to Basics: TecIO <https://www.youtube.com/watch?v=CNHONZrpeYU/>`_


ZoneType is FEPOLYGON or FEPOLYHEDRON

.. figure:: images/Polygon.png
   :width: 800
   :align: center
   
   2d case

.. figure:: images/Polygon1.png
   :width: 800
   :align: center

   node cell face figure

::
   
  Face-1  (1,2)
  Face-2  (2,3)
  Face-3  (3,4)
  Face-4  (4,5)
  Face-5  (5,6)
  Face-6  (6,1)
  Face-7  (3,7)
  Face-8  (7,8)
  Face-9  (8,9)
  Face-10 (9,10)
  Face-11 (10,4)
  Face-12 (10,11)
  Face-13 (11,12)
  Face-14 (12,13)
  Face-15 (13,5 )

  Face-1  :(1 ,2 ) Leftelem : 0 RightElem 1
  Face-2  :(2 ,3 ) Leftelem : 0 RightElem 1
  Face-3  :(3 ,4 ) Leftelem : 2 RightElem 1
  Face-4  :(4 ,5 ) Leftelem : 3 RightElem 1
  Face-5  :(5 ,6 ) Leftelem : 0 RightElem 1
  Face-6  :(6 ,1 ) Leftelem : 0 RightElem 1
  Face-7  :(3 ,7 ) Leftelem : 0 RightElem 2
  Face-8  :(7 ,8 ) Leftelem : 0 RightElem 2
  Face-9  :(8 ,9 ) Leftelem : 0 RightElem 2
  Face-10 :(9 ,10) Leftelem : 0 RightElem 2
  Face-11 :(10,4 ) Leftelem : 3 RightElem 2
  Face-12 :(10,11) Leftelem : 0 RightElem 3
  Face-13 :(11,12) Leftelem : 0 RightElem 3
  Face-14 :(12,13) Leftelem : 0 RightElem 3
  Face-15 :(13,5 ) Leftelem : 0 RightElem 3
  
.. figure:: images/Polygon2.png
   :width: 800
   :align: center

   Hexagon and 1 Octagon
   
.. figure:: images/Polygon3.png
   :width: 800
   :align: center

   Hexagon and 1 Octagon Node Face Cell Figure
     









