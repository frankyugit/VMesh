#!/usr/local/bin/python
"""VMesh Classes

This module comprises all basic classes for voxel-based mesh generation
The following example show two typical use cases of this module.

Todo:
    * Cut-Cell method
    * Octree based voexlization
    * boundary defination

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""


import datetime
from math import sqrt, hypot, acos, fabs
import struct
from gftemplates import *
import vtk



#*********************************************************
# class Point starts here
#
class Point(object):
    """The class Point defines 3D points

    Attributes:
        x (float): x coordinate
        y (float): y coordinate
        z (float): z coordinate

    The Point class maybe used as Vector class.
    Potentially Vector = Point could be added

    """

    def __init__(self, x=0, y=0, z=0):
        """Standard constructor with default value 0 for all coordinates

        Example:
            * Point() will generate (0,0,0)
            * Point(3.14) will generate (3.14,0,0)

        """
        self.x = x
        self.y = y
        self.z = z


    def __str__(self):
        """generates a string representation of a Point object
        """
        return '(%g, %g, %g)' % (self.x, self.y, self.z)

    def __add__(self, other):
        """polymorphy for the operator +

        With this method defined, we can use ''+'' to add two point objects
        p1 and p2 like in a vector addition, i.e. p1 + p2 (which is equivalent
        to p1.__add__(p2)) returns a new Point with the x coordinate determined
        by the sum of the x coordinates of self and other, and the same applies
        for the y coordinates

        Args:
            other (Point): The object on the right side of the ''+''

        Returns:
            Point: a new Point with the x coordinate determined
            by the sum of the x coordinates of p1 and p2, and the same applies
            for the y and z coordinates

        .. _See:
                http://docs.python.org/ref/specialnames.html

        """
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        """polymorphy for the operator -

        Uses the ''+'' definition above

        Args:
            other (Point): The object on the right side of the ''+''

        Returns:
            Point: a new Point with the x coordinate determined
            by the sum of the x coordinates of self (left to ''-'' object) and
            other (right to ''-'' object), and the same applies
            for the y and z coordinates

        .. _See:
                http://docs.python.org/ref/specialnames.html

        """
        return self.__add__(-1*other)

    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)


    def __mul__(self, other):
        """returns a new Point determined by multiplying coordinates independantly

        This method provide a polymophy for the ''*'' operator.
        The product self*other generates a new Point with each coordinate
        determined as product of corresponding coordinates of the calling
        object (left, self) and the point right to the ''*'', the argument
        ''other''.
        If other is not a Point a scalar is assumed and each
        coordinate of the calling object is multiplied by the scalar ''other''.

        Args:
            other (Point): The object on the right side of the ''*''

        Returns:
            Point: a new Point with the x coordinate determined
            by the product of the x coordinates of self and other, and the
            same applies for the y and z coordinates

        .. _See:
                http://docs.python.org/ref/specialnames.html

        """
        result = Point()
        try:
            result = Point(other.x * self.x, other.y * self.y, other.z * self.z)
        except (AttributeError):
            result = Point(other * self.x, other * self.y, other * self.z)
        return result

    __rmul__ = __mul__
    """implements reverse multiplication, i.e. of Point ''*''  scalar
    """

    def __div__(self, other):
        """implements division of a Point by a scalar, i.e. Point ''/'' scalar

        Returns:
            Point: each coordinate value multiplied with the inverse of the
            scalar ''other''

        """
        return self.__mul__(1.0/other)



    # With this method defined, two point objects can be compared with
    # >, <, and ==.
    def __cmp__(self, other):
        """compares 2 Points and returns -1 for self < other, 1 for self> other
        and 0 for self == other
        """
        # compare them using the x values first
        if self.x > other.x: return 1
        if self.x < other.x: return -1

        # x values are the same... check y values
        if self.y > other.y: return 1
        if self.y < other.y: return -1

        # x values are the same... check y values
        if self.z > other.z: return 1
        if self.z < other.z: return -1

        # y values are the same too. . . it's a tie
        return 0

    def __eq__(self,other):
        if self.x != other.x: return False
        if self.y != other.y: return False
        if self.z != other.z: return False
        return True

    def __hash__(self):
        return hash((self.x, self.y, self.z))


    # Other general methods
    def getCoordinates(self):
        """ returns a tuple with x, y and z coordinate
        """
        return self.x,self.y,self.z

    def getCoordinate(self,index=0):
        """ implements the index access to coordinates with x-coordinate as default
        """
        if index == 1: return self.y
        if index == 2: return self.z
        return self.x

    def setCoordinate(self,index=0, val=0.0):
        """ implements setter of the index access to coordinates with x-coordinate as default
        """
        if index == 1: self.y = val
        elif index == 2: self.z = val
        else: self.x = val


    def abs(self):
        """ returns a Point with the coordinates as abs() of calling point coodinates
        """
        return Point(abs(self.x), abs(self.y), abs(self.z))



    def distanceToLine(self,line):
        """returns distance of Point to line in argument

        Math bachground:
            The 3 points defined by Point ''self'', the orthogonal foot Point
            from Point ''self'' on the Line ''line'' and an
            auxilliary Point ''a'' on Line ''line'' form a right triangle.
            Simple trigonometry yields:
            sin(alpha) = result/hypothenuses (1) with hypothenuses=len(Line(self,a))
            Definition of scalar product:
            cos(alpha) = scalarProduct / len(Line(self,a)) / len(line) (2)
            Another trigonometric relation:
            sin(alpha) = sqrt(1 - cos(alpha)^2) (3)
            With (3) in (1) and (2) in the resulting expression for result:
            result = len(Line(self,a))*sqrt(1 - (scalarProduct / len(Line(self,a)) / len(line))^2)
            result = sqrt(len(Line(self,a))^2 - (scalarProduct / len(line))^2)
            with definitions da=len(Line(self,a)) and tmp=scalarProduct/len(line)
            result = sqrt(da^2 - tmp^2) = sqrt((da-tmp)*(da+tmp))

        """
        # Determine distances of self to line points
        d1, d2 = self.distance(line.p1), self.distance(line.p2)
        #print "d1={}, d2= {}, length p1 to p2 {}".format(d1, d2, length)
        # Determine auxiliary point a on line and distance da of self to a
        a  = line.p1+d1/(d1+d2)*(line.p2-line.p1)
        da = self.distance(a)
        if da == 0.0: return 0.0
        tmp = line.scalarProduct(Line(self,a))/line.length()
        if (da-tmp)*(da+tmp) < 1.0e-10: return 0.0
        return sqrt((da-tmp)*(da+tmp))

    def direction_from_origin(self):
        """ returns Point self scaled to unit length
        """
        if (self == origin): raise()
        return self / self.distance_from_origin()

    def distance(self, other):
        """ returns the distance of calling Point ''self'' to Point ''other''
        """
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return sqrt(dx * dx + dy * dy + dz * dz)

    def distance_from_origin(self):
        """ returns the distance of calling Point ''self'' to the origin (0,0,0)
        """
        #return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        return self.distance(origin)

    def scalarProduct(self, other):
        """ returns the scalar product ''self'' . ''other''
        """
        return (self.x * other.x + self.y * other.y + self.z*other.z)

    def vectorProduct(self, other):
        """ returns the vector product ''self'' x ''other''
        """
        return Point(   self.y * other.z - self.z * other.y,
                        self.z * other.x - self.x * other.z,
                        self.x * other.y - self.y * other.x)


    def isIn1stOctant(self):
        """ returns true if all coordinates positive
        """
        return (self.x > 0) and (self.y > 0) and (self.z > 0 )

    def getOctantPoint(self):
        """ returns the 3d octant as a unit Point, where the Point self is located in
        """
        result = Point(0,0,0)
        if   self.x >= 0 : result.x =  1
        else             : result.x = -1
        if   self.y >= 0 : result.y =  1
        else             : result.y = -1
        if   self.z >= 0 : result.z =  1
        else             : result.z = -1
        return result


    def getOctantPoint_old(self):
        result = Point(0,0,0)
        if   self.x > 0 : result.x =  1
        elif self.x < 0 : result.x = -1
        if   self.y > 0 : result.y =  1
        elif self.y < 0 : result.y = -1
        if   self.z > 0 : result.z =  1
        elif self.z < 0 : result.z = -1
        return result



    def asVertex(self, radius=0.1, color="Black"):
        return Vertex(self.x,self.y,self.z,rad,col)

    # class Point ends here
    #*********************************************************

origin = Point(0,0,0)



class Vertex(Point):
    def __init__(self, x=0, y=0, z=0, radius=0.1, color="Black"):
        """Standard constructor for Vertices which are Points with dimension
        radius and a color, so actually colored spheres in 3d space

        Example:
            * Vertex() will generate (0,0,0)
            * Vertex(3.14,0,0,0.1,(0 0 0)) will generate a black sphere
              with a radius 0.1 at point coordinates (3.14,0,0)

        """
        super(Vertex, self).__init__(x,y,z)
        self.radius = radius
        self.color = color



#*********************************************************
# class Line starts here
#
class Line(object):
    """The class Line defines 3D points

    Attributes:
        p1 (Point): first point
        p2 (Point): second point

    """

    def __init__(self, p1, p2):
        """Standard constructor

        Example:
            * Line(Point.origin,Point(3,4,5)) will generate a line
            starting at the origin and pointing to Point 3,4,5

        """
        self.p1=p1; self.p2=p2

    # With this method defined, two point objects can be compared with
    # >, <, and ==.
    def __cmp__(self, other):
        """compares 2 Lines by comparing corresponding Points

        Returns:
            -1 for self.p1 < other.p1,
            1 for self.p1 > other.p1,
            if self.p1 == other.p1 comparison continues for p2s
            0 for l1 == l2

        """
        # compare them using the x values first
        if self.p1 < other.p1: return -1
        # x values are the same... check y values
        if self.p1 > other.p1: return 1
        # self.p1 == other.p1
        # compare them using the x values first
        if self.p2 < other.p2: return -1
        # x values are the same... check y values
        if self.p2 > other.p1: return 1
        # lines are identical
        return 0

    def __str__(self):
        """defines a string representation of a Line object
        """
        return 'Line from {} to {}'.format(self.p1, self.p2,)

    def dx(self):
        """returns x2-x1
        """
        return self.p2.x - self.p1.x

    def dy(self):
        """returns y2-y1
        """
        return self.p2.y - self.p1.y

    def dz(self):
        """returns z2-z1
        """
        return self.p2.z - self.p1.z

    def length(self):
        """returns length of line
        """
        return (self.p1.distance(self.p2))

    def direction(self):
        """returns a line from origin in the direction of line self with unit length
        """
        return Line(Point(0,0,0),(self.p2-self.p1)/self.length())

    def getMidPoint(self):
        """returns the middle point of line self
        """
        return 0.5*(self.p1+self.p2)

    def orientation(self):
        """used for setting perspectives
        obsolete, delete after check
        """
        result  = [(0,0,0),0,0,0]
        # midpoint of edge between vert0 and vert1
        midpoint = self.getMidPoint()
        result[0]  = midpoint.getCoordinates()
        result[1] = self.length()
        # tip = translated endpoint (w/ midpoint at 0,0,0)
        tip = self.p2-midpoint
        # length of cylinder
        hyp = 0.5*result[1]
        # distance from y axis
        xzhyp = hypot(tip.x,tip.z)
        # rotate zx plane around north/south (y) axis to
        if xzhyp!=0: result[2] = acos(tip.z/xzhyp)
        if tip.x<=0: result[2] = -result[2]
        # rotate about x axis, tilting y axis to tip
        if hyp!=0: result[3] = acos(tip.y/hyp)
        return result


    def scalarProduct(self, other):
        """returns scalar product of lines self and other

        Calls the corresponding Point method
        """
        return (self.p2-self.p1).scalarProduct(other.p2-other.p1)

    def distanceToPoint(self,p):
        """returns the distance to a point p

        Calls the corresponding Point method
        """
        return p.distanceToLine(self)

    def isParallelShiftedTo(self, other):
        """returns True if self and other are same vector
        """
        return  (self.p1-self.p2 == other.p1-other.p2)


class Edge(Line):
    def __init__(self, p1, p2, radius=0.1, color=None):
        """Standard constructor for Edges which are Lines with dimension
        radius and a color, so actually colored cylinders in 3d space
        """
        super(Edge, self).__init__(p1,p2)
        self.radius = radius
        self.color = color


#*********************************************************
# class Triangle starts here
#
class Triangle(object):
    """defines the geometry primitive triangle defines by 3 Points

    Attributes:
        p1 (Point) : first Point
        p2 (Point) : second Point
        p3 (Point) : third Point

    Triangles are also called Faces in typical 3d tesselation files.
    So Face = Triangle could be added
    """

    def __init__(self, p1, p2, p3):
        """ a Triangle is constructed by 3 points
        """
        self.p1=p1; self.p2=p2; self.p3=p3

    def __eq__(self, other):
        """ compare self and other by checking equality of all three points

        With this method defined, two Triangle objects can be compared with
         >, <, and ==.

        Remark:
            Could be extended by rotational shifts of other points
            other.p1,other.p2,other.p3 = other.p3,other.p1,other.p2 (2 times)

        """
        return (self.p1 == other.p1) and (self.p2 == other.p2) and (self.p3 == other.p3)

    def __str__(self):
        """defines a string representation of a Triangle object
        """
        return 'Triangle ({}, {}, {})'.format(self.p1, self.p2, self.p3)

    def getAreaVector(self):
        """returns the orthoganal vector, Point respectively with the length
        of the area of the triangle calculated via the vector product
        """
        return 0.5*((self.p2-self.p1).vectorProduct(self.p3-self.p1))

    def getArea(self):
        """returns the area, calculated with getAreaVector() method
        """
        return self.getAreaVector().distance_from_origin()


    def direction(self):
        """returns orientation of the triangle calculated via the vector product
        """
        return ((p2-p1).vectorProduct(p3-p1)).direction_from_origin()

    def distanceToPoint(self,p):
        """to be implemented
        """
        pass

    def getPointsSorted(self):
        """ returns a list of tuples with edge length and opposite point
        with smallest length edge on position 0 and largest on position 2
        """
        return sorted( [(self.p3.distance(self.p2),self.p1),
                        (self.p1.distance(self.p3),self.p2),
                        (self.p2.distance(self.p1),self.p3)] )

    def getCenterPoint(self):
        return Point( (self.p1.x+self.p2.x+self.p3.x)/3.0,
                      (self.p1.y+self.p2.y+self.p3.y)/3.0,
                      (self.p1.z+self.p2.z+self.p3.z)/3.0 )

    def getExtendInDirection(self, xyz = 0):
        """returns a tuple containing actual extend of Triangle self and the
        line index of longest line in direction x
        """
        l1 = Line(self.p2,self.p3)
        l2 = Line(self.p3,self.p1)
        l3 = Line(self.p1,self.p2)
        if   xyz == 1 : dx = (abs(l1.dy()),abs(l2.dy()),abs(l3.dy()))
        elif xyz == 2 : dx = (abs(l1.dz()),abs(l2.dz()),abs(l3.dz()))
        else          : dx = (abs(l1.dx()),abs(l2.dx()),abs(l3.dx()))
        return max((val, idx) for (idx, val) in enumerate(dx))

    def getExtend(self):
        return (self.getExtendInDirection(0),self.getExtendInDirection(1),self.getExtendInDirection(2))

    def getExtendFast(self):
        """returns a tuple with the extend of the Triangle in x-,y- and z-direction
        """
        lx = [self.p1.x,self.p2.x,self.p3.x]
        ly = [self.p1.y,self.p2.y,self.p3.y]
        lz = [self.p1.z,self.p2.z,self.p3.z]
        return (max(lx)-min(lx),max(ly)-min(ly),max(lz)-min(lz))

    def getMaxExtend(self):
        return max(self.getExtendFast())

    def getMaxxyz(self, xyz):
        if xyz ==0 :
            return max(self.p1.x, self.p2.x, self.p3.x)
        elif xyz == 1 :
            return max(self.p1.y, self.p2.y, self.p3.y)
        elif xyz == 2 :
            return max(self.p1.z, self.p2.z, self.p3.z)
        else :
            print("please type in 1, 2 or 3 for x, y, z\n")

    def getMinxyz(self, xyz):
        if xyz == 0 :
            return min(self.p1.x, self.p2.x, self.p3.x)
        elif xyz == 1 :
            return min(self.p1.y, self.p2.y, self.p3.y)
        elif xyz == 2 :
            return min(self.p1.z, self.p2.z, self.p3.z)
        else :
            print("please type in 1, 2 or 3 for x, y, z\n")

    def checkxyinside(self, px = 0.0000001, py = 0.0000001):
        """project the 3D trangle onto 2d plane, check if a given 2d point is in it or not
        """
        if self.checkxyinsidebyoneedge(self.p1, self.p2,self.p3, px, py) :
            if self.checkxyinsidebyoneedge(self.p2, self.p3,self.p1, px, py) :
                if self.checkxyinsidebyoneedge(self.p3, self.p1,self.p2, px, py) :
                    return True
                else :
                    return False
            else :
                return False
        else :
            return False

    def checkxyinsidebyoneedge(self, tpl, tpr, tp, px = 0.0000001, py = 0.0000001):
        """sub-method used for method "checkxyinside"
            used two points of the triangle tpl, tpr to generate a line, then check third point tp
             and given point p are on the same side or not
        """
        if abs(tpr.x - tpl.x) < 1e-12:
            if abs(tpr.y - tpl.y) < 1e-12:
                return False
            else:
                if abs(tp.x - tpl.x) < 1e-12:
                    return False
                elif (tp.x - tpl.x)*(px - tpl.x) < 0:
                    return False
                else:
                    return True
        elif abs(tpr.y - tpl.y) < 1e-12:
            if abs(tp.y - tpl.y) < 1e-12:
                return False
            elif (tp.y - tpl.y)*(py - tpl.y) < 0:
                return False
            else:
                return True
        else:
            Ypre = tpr.y - (tpr.y - tpl.y)*(tpr.x - tp.x)/(tpr.x - tpl.x)
            YVpre = tpr.y - (tpr.y - tpl.y)*(tpr.x - px)/(tpr.x - tpl.x)
            if (Ypre - tp.y)*(YVpre - py) < 0 :
                return False
            else:
                return True

    def getzonplane(self, px = 0.0000001, py = 0.0000001):
        """set up the 3D plane equation of the trangle and get the intersected z value with a straight line perpendicular to x-y plane
        """
        planeA = self.p1.y*(self.p2.z - self.p3.z) + self.p2.y*(self.p3.z -self.p1.z) + self.p3.y*(self.p1.z - self.p2.z)
        planeB = self.p1.z*(self.p2.x - self.p3.x) + self.p2.z*(self.p3.x -self.p1.x) + self.p3.z*(self.p1.x - self.p2.x)
        planeC = self.p1.x*(self.p2.y - self.p3.y) + self.p2.x*(self.p3.y -self.p1.y) + self.p3.x*(self.p1.y - self.p2.y)
        planeD = -self.p1.x*(self.p2.y*self.p3.z - self.p3.y*self.p2.z) - self.p2.x*(self.p3.y*self.p1.z-self.p1.y*self.p3.z) - self.p3.x*(self.p1.y*self.p2.z - self.p2.y*self.p1.z)
        if abs(planeC)< 1e-14 :
            return float("inf")
        else :
            return ((-planeD -planeA*px - planeB*py)/planeC)

    def isParallelShiftedTo(self, other):
        """returns True if Triangle self and other are congruent but shifted

        Could be relaxed by rotations of other points
        """
        l1  = Line(self.p1,self.p2)
        lo1 = Line(other.p1,other.p2)
        l2  = Line(self.p1,self.p3)
        lo2 = Line(other.p1,other.p3)
        return ( l1.isParallelShiftedTo(lo1) and l2.isParallelShiftedTo(lo2) )

    def half(self, direction):
        """returns a list of 2 triangles which are generated by deviding
        the Triangle self in direction x,y or z
        """
        cp1 = self.p1.getCoordinates()[direction]
        cp2 = self.p2.getCoordinates()[direction]
        cp3 = self.p3.getCoordinates()[direction]
        lot = sorted([(cp1,1),(cp2,2),(cp3,3)])
        # node in the mmiddle
        node = lot[1][1]
        if node == 1:
            p23 = (Line(self.p2, self.p3)).getMidPoint()
            t1 = Triangle(self.p1,self.p2,p23)
            t2 = Triangle(self.p1,p23,self.p3)
        if node == 2:
            p31 = (Line(self.p3, self.p1)).getMidPoint()
            t1 = Triangle(self.p2,self.p3,p31)
            t2 = Triangle(self.p2,p31,self.p1)
        if node == 3:
            p12 = (Line(self.p1, self.p2)).getMidPoint()
            t1 = Triangle(self.p3,self.p1,p12)
            t2 = Triangle(self.p3,p12,self.p2)
        return [t1,t2]

    def split(self):
        """returns a list of 4 triangles which are generated with the nodes
        of the self Triangle and with the middle points of each side
        """
        p12 = (Line(self.p1, self.p2)).getMidPoint()
        p23 = (Line(self.p2, self.p3)).getMidPoint()
        p31 = (Line(self.p3, self.p1)).getMidPoint()
        t1 = Triangle(self.p1,p12,p31)
        t2 = Triangle(self.p2,p23,p12)
        t3 = Triangle(self.p3,p31,p23)
        t4 = Triangle(p12,p23,p31)
        return [t1,t2,t3,t4]


class Geometry(object):
    def __init__(self, cadfilename, scale):
        self.vertices=[]
        self.normals=[]
        self.faces=[]
        self.cadfilename = cadfilename
        self.scaling = scale

    def getVertices(self):
        return self.vertices

    def getEdges(self):
        """
        implemented in future
        :return: 
        """
        return [[]]

    def getFaces(self):
        """ returns a list of 3 integer lists, representing the respective vertix indices
        """
        result = []
        for f in self.faces: result.append(f[0])
        return result

    def getTriangles(self):
        """ returns a list of Triangles, generated by the faces. The triangles
        have real coordinates

        ToDo: replace faces in the parser directly by Triangles
        """
        result = []
        for f in self.faces:
            index = f[0]
            result.append(Triangle( self.vertices[index[0]],
                                    self.vertices[index[1]],
                                    self.vertices[index[2]] ))
        return result

    def getTotalSurface(self):
        result = 0.0
        tlist = self.getTriangles()
        for t in tlist: result += t.getArea()
        return result



    def printMetrics(self):
        print('---------------------------------------------------------------')
        print('Nr of Vertices   = {}'.format(len(self.vertices)))
        print('Nr of Faces      = {}'.format(len(self.faces)))
        print(' with total area = {}'.format(self.getTotalSurface()))


    def parseOBJFile(self, filename, swapyz=False):
        """
        Loads a Wavefront OBJ file.
         Not tested, will modify in future
        """
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []

        material = None
        for line in open(filename, "r"):
            if line.startswith('#'): continue
            values = line.split()
            if not values: continue
            if values[0] == 'v':
                v = map(float, values[1:4])
                if swapyz:
                    v = v[0], v[2], v[1]
                self.vertices.append(Vertex(v[0],v[1],v[2]))
            elif values[0] == 'vn':
                v = map(float, values[1:3])
                if swapyz:
                    v = v[0], v[2], v[1]
                self.normals.append(v)
            elif values[0] == 'vt':
                self.texcoords.append(map(float, values[1:3]))
            elif values[0] in ('usemtl', 'usemat'):
                material = values[1]
            elif values[0] == 'mtllib':
                self.mtl = MTL(values[1])
            elif values[0] == 'f':
                face = []
                texcoords = []
                norms = []
                for v in values[1:]:
                    w = v.split('/')
                    face.append(int(w[0])-1)
                    if len(w) >= 2 and len(w[1]) > 0:
                        texcoords.append(int(w[1]))
                    else:
                        texcoords.append(0)
                    if len(w) >= 3 and len(w[2]) > 0:
                        norms.append(int(w[2]))
                    else:
                        norms.append(0)
                self.faces.append((face, norms, texcoords, material))

    def parseASCSTLFile(self, filename, swapyz=False):
        """Loads a ASCII STL file. """
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []

        material = None
        stlfile = open(filename, "r")
        for line in stlfile :
            if line.startswith('solid'): continue
            values = line.split()
            if not values: continue
            if values[0] == 'outer': continue
            if values[0] == 'endloop': continue
            if values[0] == 'endfacet': continue
            if values[0] == 'vertex':
                v = list(map(float, values[1:4]))
                if swapyz:
                    v = v[0], v[2], v[1]
                self.vertices.append(Vertex(v[0],v[1],v[2]))
            elif values[0] == 'facet':
                v = list(map(float, values[2:5]))
                if swapyz:
                    v = v[0], v[2], v[1]
                self.normals.append(v)
        stlfile.close()
        assert (len(self.vertices)) % 3 == 0,'vertices not stored by 3-element array, please check\n'
        Ntri = int(len(self.vertices)/3)
        for gindex in range(Ntri) :    
            face = []
            face.append(gindex*3)
            face.append(gindex*3 + 1)
            face.append(gindex*3 + 2)
            material = 0
            self.faces.append((face, material))

    def parsebinarySTLFile(self, swapyz=False):
        """Loads a ASCII STL file. """
        self.vertices = []
        self.normals = []
        self.faces = []

        material = None
        stlfile = open(self.cadfilename, "rb")
        stltitle = stlfile.read(80)
        stltitletext = stltitle.decode('utf-8')
        print(stltitle)
        print("\n")
        trisum = stlfile.read(4)
        trisumnum = int.from_bytes(trisum, byteorder='little')
        print(trisumnum)
        print("\n")
        count = 1
        while stlfile :
            #print("facet\n")
            (v1, v2, v3) = struct.unpack("fff", stlfile.read(4 + 4 + 4))
            #print("%d %d %d" % (v1, v2, v3))
            v1 = self.scaling*v1
            v2 = self.scaling*v2
            v3 = self.scaling*v3
            if swapyz:
                self.normals.append(Vertex(v1, v3, v2))
            else:
                self.normals.append(Vertex(v1, v2, v3))
            #print("point1\n")
            (v1, v2, v3) = struct.unpack("fff", stlfile.read(4 + 4 + 4))
            v1 = self.scaling*v1
            v2 = self.scaling*v2
            v3 = self.scaling*v3
            #print("%d %d %d" % (v1, v2, v3))
            if swapyz:
                self.vertices.append(Vertex(v1, v3, v2))
            else:
                self.vertices.append(Vertex(v1, v2, v3))
            #print("point2\n")
            (v1, v2, v3) = struct.unpack("fff", stlfile.read(4 + 4 + 4))
            v1 = self.scaling*v1
            v2 = self.scaling*v2
            v3 = self.scaling*v3
            #print("%d %d %d" % (v1, v2, v3))
            if swapyz:
                self.vertices.append(Vertex(v1, v3, v2))
            else:
                self.vertices.append(Vertex(v1, v2, v3))
            #print("point3\n")
            (v1, v2, v3) = struct.unpack("fff", stlfile.read(4 + 4 + 4))
            v1 = self.scaling*v1
            v2 = self.scaling*v2
            v3 = self.scaling*v3
            #print("%d %d %d" % (v1, v2, v3))
            if swapyz:
                self.vertices.append(Vertex(v1, v3, v2))
            else:
                self.vertices.append(Vertex(v1, v2, v3))
            stlfile.read(2)
            count = count + 1
            if count > trisumnum :
                break
        stlfile.close()

        assert (len(self.vertices)) % 3 == 0,'vertices not stored by 3-element array, please check\n'
        Ntri = int(len(self.vertices)/3)
        for gindex in range(Ntri) :
            face = []
            face.append(gindex*3)
            face.append(gindex*3 + 1)
            face.append(gindex*3 + 2)
            material = 0
            self.faces.append((face, material))

    def getgeodim(self):
        """ returns a dictionary of xmin, xmax, ymin, ymax, zmin, zmax
        """
        geodim = {}
        vx = []
        vy = []
        vz = []
        for v in self.vertices :
            vx.append(v.x)
            vy.append(v.y)
            vz.append(v.z)
        geodim['xmin'] = min(vx)
        geodim['xmax'] = max(vx)
        geodim['ymin'] = min(vy)
        geodim['ymax'] = max(vy)
        geodim['zmin'] = min(vz)
        geodim['zmax'] = max(vz)
        return geodim
        

class Mesh(object):
    """The class Mesh defines a computational mesh consisting of the
    following attributes:
        coords : list of x-, y- and z-coordinates
        Voxelised : 1d array of the filled condition for every cell 


    A Mesh maybe written, stored respectively, in the OF or NetCDF format.

    """
    def __init__(self, nodecoords, geo):
        """ standard constructor with per default empty coordinates.

        As Mesh is a visual object by inheritance the constructor calls the constructor
        of the visual qualities, sets default colours and display properties
        regarding nodes, lines and faces.

        TODO: The default values should be changed to at least containing
        2 values, e.g. xcoords=[0,1,2,3,4,5,6,7,8,9],ycoords[0,1],zcoords=[0,1]
        """
        self.nodecoords = nodecoords
        self.voxel3dwidthlist = [[],[],[]]
        self.set3DVoxelWidthList()
        self.geo = geo
        self.Voxelised = []

    def __str__(self):
       return '(x: {} \n y: {}; \n z: {})'.format(
                ','.join(str(e) for e in self.coords[0]),
                ','.join(str(e) for e in self.coords[1]),
                ','.join(str(e) for e in self.coords[2]) )

    def set3DVoxelWidthList(self):
        self.voxel3dwidthlist = [[],[],[]]
        for dimxyz in range(3):
            for Nodei in range(1, self.getNNodesxyz(dimxyz)):
                self.voxel3dwidthlist[dimxyz].append(self.nodecoords[dimxyz][Nodei] - self.nodecoords[dimxyz][Nodei-1])

    def getVoxelCenterCoords(self):
        voxelcentercoords = [[],[],[]]
        for dimxyz in range(3):
            for Voxeli in range(self.getNVoxelxyz(dimxyz)):
                voxelcentercoords[dimxyz].append(self.nodecoords[dimxyz][Voxeli] + 0.5*self.voxel3dwidthlist[dimxyz][Voxeli])
        return voxelcentercoords


    def getNNodesxyz(self, xyz):
        return len(self.nodecoords[xyz])

    def getNVoxelxyz(self, xyz):
        return (len(self.nodecoords[xyz]) - 1)

    def getNodeCoordx(self,ix):
        return self.nodecoords[0][ix]

    def getNodeCoordy(self,iy):
        return self.nodecoords[1][iy]

    def getNodeCoordz(self,iz):
        return self.nodecoords[2][iz]

    def getNodexyzCoords(self, xyz):
        return self.nodecoords[xyz]

    def getVoxelxyzWidthList(self, xyz):
        return self.voxel3dwidthlist[xyz]

    def getVoxelWidthx(self,ix):
        return self.getVoxelxyzWidthList(0)[ix]

    def getVoxelWidthy(self,iy):
        return self.getVoxelxyzWidthList(1)[iy]

    def getVoxelWidthz(self,iz):
        return self.getVoxelxyzWidthList(2)[iz]

    def get1dNodeCoords(self):
        v = []
        for x in self.nodecoords[2]:
            for y in self.nodecoords[1]:
                for z in self.nodecoords[0]:
                    v.append(Vertex(x,y,z))
        return v

    def getMeshdim(self):
        """ returns a dictionary of xmin, xmax, ymin, ymax, zmin, zmax
        """
        meshdim = {}
        meshdim['xmin'] = self.getNodeCoordx(0)
        meshdim['xmax'] = self.getNodeCoordx(-1)
        meshdim['ymin'] = self.getNodeCoordy(0)
        meshdim['ymax'] = self.getNodeCoordy(-1)
        meshdim['zmin'] = self.getNodeCoordz(0)
        meshdim['zmax'] = self.getNodeCoordz(-1)
        return meshdim

    def checkMesh(self):
        """
        check if any coord of x, y, or z exceeds meshdim
        check if all of the three arrays, x, y and z, are increasing from low to high index
        :return: 
        """
        print("Mesh Checking\n")
        print("---------------------------------------------------------------")
        meshdim = self.getMeshdim()
        Xmin = meshdim['xmin']
        Xmax = meshdim['xmax']
        Ymin = meshdim['ymin']
        Ymax = meshdim['ymax']
        Zmin = meshdim['zmin']
        Zmax = meshdim['zmax']
        XCoord = self.getNodexyzCoords(0)
        YCoord = self.getNodexyzCoords(1)
        ZCoord = self.getNodexyzCoords(2)
        print("Coordinate checking")
        for xi in range(1,len(XCoord)):
            if(XCoord[xi] < XCoord[xi-1]):
                print("Fatal error: x coordinate %d is large than %d" % (XCoord[xi-1], XCoord[xi]))
        for yi in range(1,len(YCoord)):
            if(YCoord[yi] < YCoord[yi-1]):
                print("Fatal error: y coordinate %d is large than %d" % (YCoord[yi-1], YCoord[yi]))
        for zi in range(1,len(ZCoord)):
            if(ZCoord[zi] < ZCoord[zi-1]):
                print("Fatal error: z coordinate %d is large than %d" % (ZCoord[zi-1], ZCoord[zi]))
        print("Dim checking\n")
        if Xmin != min(XCoord):
            print("Fatal error: x coordinate exceeds the minimum of mesh dim")
        if Xmax != max(XCoord):
            print("Fatal error: x coordinate exceeds the maximum of mesh dim")
        if Ymin != min(YCoord):
            print("Fatal error: y coordinate exceeds the minimum of mesh dim")
        if Ymax != max(YCoord):
            print("Fatal error: y coordinate exceeds the maximum of mesh dim")
        if Zmin != min(ZCoord):
            print("Fatal error: z coordinate exceeds the minimum of mesh dim")
        if Zmax != max(ZCoord):
            print("Fatal error: z coordinate exceeds the maximum of mesh dim")

    def setNodeXYZCoordsToList(self, xyz, nodexyzcoords):
        """sets coordinates in direction xyz """
        self.nodecoords[xyz]= nodexyzcoords

    def addNodexyzCoord(self, xyz, c):
        self.nodecoords[xyz].append(c)
        self.nodecoords[xyz].sort()

    def addNodeCoordOfPoint(self, p):
        self.addNodexyzCoord(0,p.x);self.addNodexyzCoord(1,p.y);self.addNodexyzCoord(2,p.z)

    def getPointwithIndex(self, i,j,k):
        return Point(self.nodecoords[0][i],self.nodecoords[1][j],self.nodecoords[2][k])

    def getPointwithIndexPoint(self, ip):
        return Point(self.nodecoords[0][ip.x],self.nodecoords[1][ip.y],self.nodecoords[2][ip.z])

    def getMinNodeCoord(self, xyz):
        return min(self.nodecoords[xyz])

    def getMaxNodeCoord(self, xyz):
        return max(self.nodecoords[xyz])

    def getNVoxel(self):
        return (self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*self.getNVoxelxyz(2))

    def getMeshField(self):
        return self.Voxelised

    def getLength(self, xyz):
        return (self.getMaxNodeCoord(xyz)-self.getMinNodeCoord(xyz))

    def getTotalVolume(self):
        return (self.getLength(0)*self.getLength(1)*self.getLength(2))

    def getMidPoint(self):
        px = self.getNodeCoordx(0)+self.getNodeCoordx(-1)
        py = self.getNodeCoordy(0)+self.getNodeCoordy(-1)
        pz = self.getNodeCoordz(0)+self.getNodeCoordz(-1)
        result = Point(px, py, pz)
        return (0.5 * result)

    def getVoxelxyzindexfromGlobal(self, voxeli):
        """
        mapping index of 1d voxel array to index of 3d array
        :param voxeli: 
        :return: 
        """
        pass


    def getVoxelGrowths(self, xyz):
        result = list()
        widths = self.voxel3dwidthlist[xyz]
        prev = widths[0]
        for i in range(1,len(widths)):
            current = widths[i]
            result.append(current/prev)
            prev = current
        return result

    def printMetrics(self):
        print(self.getMetricsString())


    def getMetricsString(self):
        result = ""
        result += 'Coordinates nx = {}, ny = {}, nz = {}'\
        .format(self.getNNodesxyz(0), self.getNNodesxyz(1),self.getNNodesxyz(2))
        result += 'Number of Cells nx = {}, ny = {}, nz = {}; n_total = {}'\
        .format(self.getNVoxelxyz(0), self.getNVoxelxyz(1),self.getNVoxelxyz(2),self.getNVoxel())
        result += 'Volume = length x depth x height = {} x {} x {} = {}'\
        .format(self.getLength(0),self.getLength(1),self.getLength(2),self.getTotalVolume())
        result += '---------------------------------------------------------------\n'
        result += 'Min and Max cell width in'
        tmp = self.getVoxelxyzWidthList(0)
        result += 'x-direction:  {} - {} '.format(min(tmp),max(tmp))
        tmp = self.getVoxelxyzWidthList(1)
        result += 'y-direction:  {} - {} '.format(min(tmp),max(tmp))
        tmp = self.getVoxelxyzWidthList(2)
        result += 'z-direction:  {} - {} '.format(min(tmp),max(tmp))
        result += '---------------------------------------------------------------\n'
        return result

    def getVoxel3DIndexFromCoord(self, coordx, coordy, coordz):
        """
        search the voxel x y z index for a given location
        Note: the returned index all starting from 1
        :param coordx: given x coordinate
        :param coordy: given y coordinate
        :param coordz: given z coordinate
        :return: voxel indices
        """
        indices = [0,0,0]
        coords = [coordx, coordy, coordz]
        for direction in range(3):
            meshcoord = self.getNodexyzCoords(direction)
            if coords[direction] < meshcoord[0]:
                indices[direction] = -1
                print("Given c coordinate exceeds the lower limit of the mesh")
            else:
                for Nodei in range(1, self.getNNodesxyz(direction)):
                    if coords[direction] < meshcoord[Nodei]:
                        indices[direction] = Nodei
                        break
                if coords[direction] > meshcoord[-1]:
                    indices[direction] = -1
                    print("Given c coordinate exceeds the upper limit of the mesh")
        return indices

    def getVoxelCenterCoordFromIndices(self, iX, iY, iZ):
        coords = [0.0,0.0,0.0]
        voxcentercoord = self.getVoxelCenterCoords()
        indexX = iX - 1
        indexY = iY - 1
        indexZ = iZ - 1
        indices = [indexX, indexY, indexZ]

        for direction in range(3):
            index = indices[direction]
            if index < 0 or index > (self.getNVoxelxyz(0) - 1):
                coords[direction] = -1
                print("Given index exceeds the range of indices")
            else:
                coords[direction] = voxcentercoord[direction][index]

        return coords

    # def getIndex(self, xyz ,c):
    #     i=0; found = False
    #     while not found and i<len(self.coords[xyz]):
    #         if self.coords[xyz][i] >= c: found = True
    #         else: i += 1
    #     if i==0: return 0
    #     if not found: return len(self.coords[xyz])-1
    #     if c-self.coords[xyz][i-1] < self.coords[xyz][i]-c : return i-1
    #     return i

    # def getIndexOfPoint(self, p):
    #     return Point(self.getIndex(0,p.x),self.getIndex(1,p.y),self.getIndex(2,p.z))


    # def isValidIndexPoint(self, ip):
    #     if ip.x < 0 : return False
    #     if ip.y < 0 : return False
    #     if ip.z < 0 : return False
    #     if ip.x >= len(self.coords[0]) : return False
    #     if ip.y >= len(self.coords[1]) : return False
    #     if ip.z >= len(self.coords[2]) : return False
    #     return True

    def setMeshFieldiValue(self,i, j, k, value):
        """
        artifically set up certain element of field to given value
        :param i: index of x starting from 0
        :param j: index of y starting from 0
        :param k: index of z starting from 0
        :param value:  given field value
        """
        index = self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*j + i
        self.Voxelised[index] = value

    def setVoxelising(self):
        """based on ray scan to each of the trangle
        get the filled condition of each voxel
        """
        start = datetime.datetime.now()
        triangleList = self.geo.getTriangles()
        for k in range(self.getNVoxelxyz(2)) :
            for j in range(self.getNVoxelxyz(1)):
                for i in range(self.getNVoxelxyz(0)):
                    self.Voxelised.append(0)

        voxelcentercoords = self.getVoxelCenterCoords()
        xcoords = voxelcentercoords[0]
        ycoords = voxelcentercoords[1]
        zcoords = voxelcentercoords[2]
        Nx = self.getNVoxelxyz(0)
        Ny = self.getNVoxelxyz(1)
        Nz = self.getNVoxelxyz(2)
        Nodezcoords = self.getNodexyzCoords(2)
        NNodez = self.getNNodesxyz(2)
        ###possible increase speed with Numpy###
        crosstList = []
        possiblecrosstList = []
        possiblecrosstListy = []
        print("total Cell = %d\n" % (Nx*Ny*Nz))
        print("xmin = %f\n" % (xcoords[0]))
        print("xmax = %f\n" % (xcoords[Nx - 1]))
        print("ymin = %f\n" % (ycoords[0]))
        print("ymax = %f\n" % (ycoords[Ny - 1]))
        print("zmin = %f\n" % (zcoords[0]))
        print("zmax = %f\n" % (zcoords[Nz - 1]))
        for j in range(Ny):
            possiblecrosstListy = []
            for t in triangleList:
                if t.getMinxyz(1) < ycoords[j] and t.getMaxxyz(1) > ycoords[j]:
                    possiblecrosstListy.append(t)
            for i in range(Nx):
                possiblecrosstList = []
                count=0
                for ty in possiblecrosstListy:
                    count = count+1
                    if ty.getMinxyz(0) < xcoords[i] and ty.getMaxxyz(0) > xcoords[i]:
                        if ty.checkxyinside(xcoords[i], ycoords[j]):
                            if ty.getzonplane(xcoords[i], ycoords[j]) > (self.getNodeCoordz(-1) + 1e-6):
                                print("Found face parllel to xy-plane at x = %d, y = %d\n" % (i, j))
                                print("Check the closeness of geometry\n")
                                continue
                            elif ty.getzonplane(xcoords[i], ycoords[j]) < (self.getNodeCoordz(0) - 1e-6):
                                print("Found face parllel to xy-plane at x = %d, y = %d\n" % (i, j))
                                print("Check the closeness of geometry\n")
                                continue
                            else:
                                possiblecrosstList.append(ty)
                if len(possiblecrosstList) == 0 :
                    continue
                elif len(possiblecrosstList) == 1 :
                    continue
#                if len(possiblecrosstList) % 2 :
#                    print "Found odd number of faces on one ray at x = %d, y = %d\n" % (i, j)
#                    print "Check the closeness of geometry\n"
#                    continue
                else :
                    Zpossible = []
                    for tp in possiblecrosstList :
                        intersectZ = tp.getzonplane(xcoords[i], ycoords[j])
                        repeatedZ = False
                        for zpossiblei in Zpossible:
                            if abs(intersectZ - zpossiblei) > 1e-14:
                                pass
                            else:
                                repeatedZ = True
                        if repeatedZ:
                            pass
                        else:
                            Zpossible.append(intersectZ)
                    Zindices = []
                    for zcoord in Zpossible:
                        if zcoord < Nodezcoords[0] or zcoord > Nodezcoords[-1]:
                            print("Found intercepted z beyond the scope of mesh\n")
                        for zi in range(NNodez):
                            if zcoord < Nodezcoords[zi]:
                                Zindices.append(zi-1)
                                break
                    Zindices = sorted(Zindices)
                    if len(Zindices)%2:
                        print("There is odd number of intecepted z coordinates on (%d, %d), please check" % (i,j))
                    pairs = len(Zindices)//2
                    for pairindex in range(pairs) :
                        Zmin = Zindices[2*pairindex]
                        Zmax = Zindices[2*pairindex + 1]
                        for Zindex in range(Zmin, Zmax):
                            self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*Zindex + self.getNVoxelxyz(0)*j + i] = 1

        # mark = 2
        # for k in range(self.getNNodes(2)) :
        #     for j in range(self.getNNodes(1)):
        #         for i in range(self.getNNodes(0)):
        #             mark = self.innercheck(i, j, k, mark)


        end = datetime.datetime.now()
        print("Total excuted time %f" % (end-start).seconds)

    def innercheck(self, i = 0, j = 0, k = 0, mark = 2, recuriter = 0):
        print("i=%d, j=%d, k=%d\n"%(i, j, k))
        mark = mark
        recursiontimes = recuriter
        if recursiontimes > 500:
            print("Recursion exceeds 500 cells, either there's no independent structure or it's too large")
            return -1
        if i > 0:
            if self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*j + i-1] == 1:
                pass
            elif self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*j + i-1] > (mark-1):
                pass
            else:
                self.Voxelised[self.getNVoxelxyz(0) * self.getNVoxelxyz(1) * k + self.getNVoxelxyz(0) * j + i-1] = mark
                recursiontimes = recursiontimes + 1
                self.innercheck(i-1, j ,k, mark, recursiontimes)
        else:
            pass
        if j > 0:
            if self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*(j-1) + i] == 1:
                pass
            elif self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*(j-1) + i] > (mark-1):
                pass
            else:
                self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*(j-1) + i] = mark
                recursiontimes = recursiontimes + 1
                self.innercheck(i, j-1 ,k, mark, recursiontimes)
        else:
            pass
        if k > 0:
            if self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-1) + self.getNVoxelxyz(0)*j + i] == 1:
                pass
            elif self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-1) + self.getNVoxelxyz(0)*j + i] > (mark-1):
                pass
            else:
                self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-1) + self.getNVoxelxyz(0)*j + i] = mark
                recursiontimes = recursiontimes + 1
                self.innercheck(i, j ,k-1, mark, recursiontimes)
        else:
            pass
        if i < (self.getNVoxelxyz(0)-1):
            if self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*j + i+1] == 1:
                pass
            elif self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*j + i+1] > (mark-1):
                pass
            else:
                self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0) * j + i+1] = mark
                recursiontimes = recursiontimes + 1
                self.innercheck(i+1, j ,k, mark, recursiontimes)
        else:
            pass
        if j < (self.getNVoxelxyz(1)-1):
            if self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*(j+1) + i] == 1:
                pass
            elif self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*(j+1) + i] > (mark-1):
                pass
            else:
                self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*k + self.getNVoxelxyz(0)*(j+1) + i] = mark
                recursiontimes = recursiontimes + 1
                self.innercheck(i, j+1 ,k, mark, recursiontimes)
        else:
            pass
        if k < (self.getNVoxelxyz(2)-1):
            if self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k+1) + self.getNVoxelxyz(0)*j + i] == 1:
                pass
            elif self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k+1) + self.getNVoxelxyz(0)*j + i] > (mark-1):
                pass
            else:
                self.Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k+1) + self.getNVoxelxyz(0)*j + i] = mark
                recursiontimes = recursiontimes + 1
                self.innercheck(i, j ,k+1, mark, recursiontimes)
        else:
            pass
            # if not recursiontimes:
            #     print("Recursion in %d, %d, %d for %d times"%(i, j, k, recursiontimes))
            #     mark = mark + 1
            # return mark

    def cleanMesh(self):
        """ to be implemented in future
        """
        pass


    def parseGFFile(self, filename):
        """parses a GASFLOW file for coordinates, walls, obstacles and holes """
        self.__init__()
        isInMeshgn = False; indXYZ = -1
        isInXput = False

        for rawline in open(filename, "r"):
            #print rawline
            # remove comments introduced by leading ';'
            line = rawline.split(';',1)[0]
            # remove inparticular leading blanks
            line = line.strip()
            # skip empty lines
            if line == "": continue

            if line.startswith('$meshgn'):
                #print "$meshgn found in file"
                isInMeshgn = True; self.coords=[[],[],[]]; continue

            if isInMeshgn:
                if line.startswith('$end'):
                    isInMeshgn = False; indXYZ = -1; continue

                cleanline = line.replace("="," ")
                cleanline = cleanline.replace(","," ")
                values = cleanline.split()
                if not values: continue

                if   values[0]=='xgrid': indXYZ=0; del values[0]
                elif values[0]=='ygrid': indXYZ=1; del values[0]
                elif values[0]=='zgrid': indXYZ=2; del values[0]

                if indXYZ > -1: self.coords[indXYZ].extend(map(float,values))

            ### xput contains obstacles, walls and holes ###

            if line.startswith('$xput'):
                # print "$xput found in file"
                isInXput = True
                self.walls = []; self.obstacles=[]; self.holes = []
                continue

            if isInXput:
                if line.startswith('$end'):
                    isInXput = False; continue

                cleanline = line.replace("=", " ")
                cleanline = cleanline.replace(",", " ")
                cleanline = cleanline.replace("(", " ")
                cleanline = cleanline.replace(")", " ")
                cleanline = cleanline.replace(":", " ")
                values = cleanline.split()
                #if not len(values)== 12: print "!!! Line with strange content found:\n{}".format(line); continue
                #   0   1 2 3    4   5  6  7    8   9  10  11
                # walls(1:8,1) = 1, 43, 1, 41, 33, 33, 1, 1,
                if values[0]=='mobs':
                    (ix1,ix2,iy1,iy2,iz1,iz2) = [int(i)-1 for i in values[4:10]]
                    self.obstacles.append(Obstacle(ix1,ix2,iy1,iy2,iz1,iz2))
                if values[0]=='walls':
                    (ix1,ix2,iy1,iy2,iz1,iz2) = [int(i)-1 for i in values[4:10]]
                    self.walls.append(Wall(ix1,ix2,iy1,iy2,iz1,iz2))
                if values[0]=='holes':
                    (ix1,ix2,iy1,iy2,iz1,iz2) = [int(i)-1 for i in values[4:10]]
                    self.holes.append(Obstacle(ix1,ix2,iy1,iy2,iz1,iz2))

    def writeOFFile(self, filename, pathname ='./'):
        Voxelised = self.Voxelised
        print("OF File Writing...")
        fid = open(pathname + filename, 'w')
        fid.write("FoamFile\n")
        fid.write("{\n")
        fid.write("\tversion\t2.0;\n")
        fid.write("\tformat\tascii;\n")
        fid.write("\tclass\tvolScalarField;\n")
        fid.write("\tlocation\t\"0\";\n")
        fid.write("\tobject\tCell;\n")
        fid.write("}\n")
        fid.write("dimensions\t[0 0 0 0 0 0 0];\n\n")
        fid.write("internalField\tnonuniform List<scalar>\n")
        fid.write("%d\n" % (self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*self.getNVoxelxyz(2)))
        fid.write("(\n")
        for k in range(self.getNVoxelxyz(2)) :
            for j in range(self.getNVoxelxyz(1)) :
                for i in range(self.getNVoxelxyz(0)) :
                    fid.write("%d\n" % (Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-1) + self.getNVoxelxyz(0)*(j-1) + i]))
        fid.write(")\n;\n\n")
        fid.write("boundaryField\n")
        fid.write("{\n")
        fid.write("\twall\n")
        fid.write("\t{\n")
        fid.write("\t\ttype\tfixedValue;\n")
        fid.write("\t\tvalue\tuniform 1;\n")
        fid.write("\t}\n")
        fid.write("}\n")
        fid.close()
        print("Done")

    def writeNetCDFFile(self, filename, pathname ='./'):
        Voxelised = self.Voxelised
        Voxelisedwithg = []
        xgrid = []
        ygrid = []
        zgrid = []
        NCellx = self.getNVoxelxyz(0) + 2
        NCelly = self.getNVoxelxyz(1) + 2
        NCellz = self.getNVoxelxyz(2) + 2
        TotalCellsnum = NCellx*NCelly*NCellz
        NNodex = self.getNNodesxyz(0) + 2
        NNodey = self.getNNodesxyz(1) + 2
        NNodez = self.getNNodesxyz(2) + 2
        for k in range(NCellz):
            for j in range(NCelly):
                for i in range(NCellx):
                    if k == 0 or k == (NCellz - 1):
                        Voxelisedwithg.append(-1)
                    elif j == 0 or j == (NCelly - 1):
                        Voxelisedwithg.append(-1)
                    elif i == 0 or i == (NCellx - 1):
                        Voxelisedwithg.append(-1)
                    else :
                        if Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-2) + self.getNVoxelxyz(0)*(j-2) + i-1] == 0:
                            Voxelisedwithg.append(1)
                        elif Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-2) + self.getNVoxelxyz(0)*(j-2) + i-1] == 1:
                            Voxelisedwithg.append(0)
                        else:
                            Voxelisedwithg.append(Voxelised[self.getNVoxelxyz(0)*self.getNVoxelxyz(1)*(k-2) + self.getNVoxelxyz(0)*(j-2) + i-1])

        print("xorign = %f" % (self.getNodeCoordx(0) - self.getVoxelWidthx(0)))
        print("yorign = %f" % (self.getNodeCoordy(0) - self.getVoxelWidthy(0)))
        print("zorign = %f" % (self.getNodeCoordz(0) - self.getVoxelWidthz(0)))
        print("NetCDF File Writing...")
        xgrid.append(self.getNodeCoordx(0) - self.getVoxelWidthx(0))
        for i in range(self.getNNodesxyz(0)):
            xgrid.append(self.getNodeCoordx(i))
        xgrid.append(self.getNodeCoordx(-1) + self.getVoxelWidthx(-1))
        ygrid.append(self.getNodeCoordy(0) - self.getVoxelWidthy(0))
        for j in range(self.getNNodesxyz(1)):
            ygrid.append(self.getNodeCoordy(j))
        ygrid.append(self.getNodeCoordy(-1) + self.getVoxelWidthy(-1))
        zgrid.append(self.getNodeCoordz(0) - self.getVoxelWidthz(0))
        for k in range(self.getNNodesxyz(2)):
            zgrid.append(self.getNodeCoordz(k))
        zgrid.append(self.getNodeCoordz(-1) + self.getVoxelWidthz(-1))
        NetCDFdict = dict(
            Totalcells = str(TotalCellsnum),
            xnodes=str(NNodex),
            ynodes=str(NNodey),
            znodes=str(NNodez),
            xcells     = str(NCellx),
            ycells     = str(NCelly),
            zcells     = str(NCellz),
            Voxelisedwithg  = str(Voxelisedwithg).replace('[','').replace(']',''),
            xgrid           = str(xgrid).replace('[','').replace(']',''),
            ygrid           = str(ygrid).replace('[','').replace(']',''),
            zgrid           = str(zgrid).replace('[','').replace(']',''))
        fid = open(pathname + filename, 'w')
        fid.write(gf_netcdftemplate.substitute(NetCDFdict))
        fid.close()
        print("Done\n")


class VTKScene (object) :
    """The class VTKScene defines a scene to render the mesh, as well as the field inside it, by VTK
       it consists of the following attributes:
        coords : list of x-, y- and z-coordinates



    """
    def __init__(self, rmesh, fieldlist, meshopacity, VTKrenderWI):
        self.xCoords = vtk.vtkFloatArray()
        self.yCoords = vtk.vtkFloatArray()
        self.zCoords = vtk.vtkFloatArray()
        self.mesh = rmesh
        self.Xmin = self.mesh.getMeshdim()['xmin']
        self.Xmax = self.mesh.getMeshdim()['xmax']
        self.Ymin = self.mesh.getMeshdim()['ymin']
        self.Ymax = self.mesh.getMeshdim()['ymax']
        self.Zmin = self.mesh.getMeshdim()['zmin']
        self.Zmax = self.mesh.getMeshdim()['zmax']

        self.meshopacity = meshopacity
        self.CellNumX = self.mesh.getNVoxelxyz(0)
        self.CellNumY = self.mesh.getNVoxelxyz(1)
        self.CellNumZ = self.mesh.getNVoxelxyz(2)
        self.field = self.mesh.getMeshField()
        self.RTGrid = vtk.vtkRectilinearGrid()
        self.FieldToView = fieldlist
        self.VTKrenWI = VTKrenderWI
        self.ren = self.VTKrenWI.GetRenderWindow().GetRenderers().GetFirstRenderer()
        print(self.VTKrenWI.GetRenderWindow().GetRenderers())

    def setVTKSceneCoords(self):
        NNodex = self.CellNumX + 1
        for i in range(NNodex):
            self.xCoords.InsertNextValue(self.mesh.getNodeCoordx(i))
        NNodey = self.CellNumY + 1
        for j in range(NNodey):
            self.yCoords.InsertNextValue(self.mesh.getNodeCoordy(j))
        NNodez = self.CellNumZ + 1
        for k in range(NNodez):
            self.zCoords.InsertNextValue(self.mesh.getNodeCoordz(k))

    def vtkRTGridGen(self):
        NNodex = self.CellNumX + 1
        NNodey = self.CellNumY + 1
        NNodez = self.CellNumZ + 1
        self.RTGrid.SetDimensions(NNodex, NNodey, NNodez)
        self.RTGrid.SetXCoordinates(self.xCoords)
        self.RTGrid.SetYCoordinates(self.yCoords)
        self.RTGrid.SetZCoordinates(self.zCoords)
        scalarfield = vtk.vtkFloatArray()
        for k in range(self.CellNumZ) :
            for j in range(self.CellNumY) :
                for i in range(self.CellNumX) :
                    scalarfield.InsertNextValue(self.field[self.CellNumX*self.CellNumY*k + self.CellNumX*j + i])
        # print(scalarfield)
        self.RTGrid.GetCellData().SetScalars(scalarfield)

    def solidVisualOn(self):
        self.SolidOrGasVisual = True

    def gasVisualOn(self):
        self.SolidOrGasVisual = False

    def setFieldToView(self, fstartindex,fendindex):
        """
        set up the fields going to be viewed by providing filed indices
        :param fstartindex: start field index
        :param fendindex: end field index
        """
        self.FieldToView["start"] = fstartindex
        self.FieldToView["end"] = fendindex

    def vtkRender(self):

        self.setVTKSceneCoords()
        self.vtkRTGridGen()

        threshold = vtk.vtkThreshold()
        threshold.SetInputData(self.RTGrid)
        threshold.ThresholdBetween(self.FieldToView["start"], self.FieldToView["end"])
        threshold.Update()

        GridMapper = vtk.vtkDataSetMapper()
        GridMapper.SetInputConnection(threshold.GetOutputPort())
        GridActor = vtk.vtkActor()
        #GridActor.VisibilityOn()
        GridActor.SetMapper(GridMapper)
        GridActor.GetProperty().SetRepresentationToSurface()
        GridActor.GetProperty().SetColor(0.2, 0.8, 0.6)
        GridActor.GetProperty().SetOpacity(self.meshopacity)
        GridActor.GetProperty().RenderPointsAsSpheresOn()
        GridActor.GetProperty().SetEdgeColor(0.6, 0.6, 0.6)
        GridActor.GetProperty().EdgeVisibilityOn()

        cubeAxesActor = vtk.vtkCubeAxesActor()
        cubeAxesActor.SetBounds(self.RTGrid.GetBounds())
        cubeAxesActor.SetCamera(self.ren.GetActiveCamera())
        cubeAxesActor.GetTitleTextProperty(0).SetColor(1.0, 0.0, 0.0)
        cubeAxesActor.GetTitleTextProperty(0).SetColor(1.0, 0.0, 0.0)
        cubeAxesActor.GetLabelTextProperty(0).GetVerticalJustification()
        cubeAxesActor.GetTitleTextProperty(0).ShadowOn()
        cubeAxesActor.GetTitleTextProperty(0).BoldOn()
        cubeAxesActor.GetTitleTextProperty(1).SetColor(0.0, 1.0, 0.0)
        cubeAxesActor.GetLabelTextProperty(1).SetColor(0.0, 1.0, 0.0)
        cubeAxesActor.GetLabelTextProperty(1).SetFontSize(8)
        cubeAxesActor.GetTitleTextProperty(2).SetColor(0.0, 0.0, 1.0)
        cubeAxesActor.GetLabelTextProperty(2).SetColor(0.0, 0.0, 1.0)
        cubeAxesActor.SetLabelScaling(True, 100, 100, 100)
        cubeAxesActor.SetAxisBaseForY(0,100,0)
        cubeAxesActor.DrawXGridlinesOn()
        cubeAxesActor.DrawYGridlinesOn()
        cubeAxesActor.DrawZGridlinesOn()
        cubeAxesActor.SetXAxisRange(self.Xmin, self.Xmax)
        cubeAxesActor.YAxisLabelVisibilityOn()
        cubeAxesActor.ZAxisMinorTickVisibilityOn()
        cubeAxesActor.SetGridLineLocation(cubeAxesActor.VTK_GRID_LINES_FURTHEST)
        #cubeAxesActor.XAxisMinorTickVisibilityOff()
        #cubeAxesActor.YAxisMinorTickVisibilityOff()
        #cubeAxesActor.ZAxisMinorTickVisibilityOff()

        # reader = vtk.vtkSTLReader()
        # reader.SetFileName("bCtl.stl")
        # reader.Update()
        # print(reader.GetFileName())
        # sMapper = vtk.vtkDataSetMapper()
        # sMapper.SetInputConnection(reader.GetOutputPort())
        # sActor = vtk.vtkLODActor()
        # sActor.SetMapper(sMapper)
        # sActor.GetProperty().SetRepresentationToSurface()
        # sActor.GetProperty().SetColor(light_grey)
        self.ren.AddActor(cubeAxesActor)
        #ren.AddActor(sActor)
        self.ren.AddActor(GridActor)
        self.ren.SetBackground(0.1, 0.2, 0.4)
        self.ren.ResetCamera()
        self.ren.GetActiveCamera().Zoom(1.5)

        self.VTKrenWI.Initialize()
        self.VTKrenWI.Render()
        self.VTKrenWI.Start()


    def renderpickupcell(self, i = 0, j = 0, k = 0):
        pickupcellid = self.RTGrid.ComputeCellId([i, j ,k])
        print("cell id = %d" % pickupcellid)
        pickupcell = self.RTGrid.GetCell(pickupcellid)
        print(pickupcell.GetPoints())
        pickupcell.GetPointIds().SetId(0, 0)
        pickupcell.GetPointIds().SetId(1, 1)
        pickupcell.GetPointIds().SetId(2, 2)
        pickupcell.GetPointIds().SetId(3, 3)
        pickupcell.GetPointIds().SetId(4, 4)
        pickupcell.GetPointIds().SetId(5, 5)
        pickupcell.GetPointIds().SetId(6, 6)
        pickupcell.GetPointIds().SetId(7, 7)
        pickupug = vtk.vtkUnstructuredGrid()
        pickupug.SetPoints(pickupcell.GetPoints())
        pickupug.InsertNextCell(pickupcell.GetCellType(), pickupcell.GetPointIds())
        pickupMapper = vtk.vtkDataSetMapper()
        pickupMapper.SetInputData(pickupug)
        pickupActor = vtk.vtkActor()
        pickupActor.SetMapper(pickupMapper)
        pickupActor.GetProperty().SetColor(0.8, 0.8, 0.0)
        pickupActor.GetProperty().SetRepresentationToSurface()
        pickupActor.GetProperty().SetOpacity(1.0)
        pickupActor.GetProperty().RenderPointsAsSpheresOn()
        pickupActor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
        pickupActor.GetProperty().SetLineWidth(3)
        pickupActor.GetProperty().EdgeVisibilityOn()
        self.ren.AddActor(pickupActor)
        self.VTKrenWI.GetRenderWindow().AddRenderer(self.ren)

    def innercheck(self, i = 0, j = 0, k = 0, mark = 2):
        self.mesh.setMeshFieldiValue(i, j, k, mark)
        self.mesh.innercheck(i,j,k,mark)
        print("Above cells have been filled by field %d" % mark)

    def vtkwrite(self):
        # write the voxelising mesh into ascii vtk file

        vtkwriter = vtk.vtkRectilinearGridWriter()
        vtkwriter.SetFileName("voxelising.vtk")
        vtkwriter.SetInputData(self.RTGrid)
        vtkwriter.SetFileTypeToASCII()
        vtkwriter.Write()












