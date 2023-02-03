import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Random;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.imageio.ImageIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import javax.swing.*;
import java.awt.*;

public class Main{

    private static final Random random = new Random();

    public static void main(String[] args) {

        int width = 1000, height = 1000;

        Screen screen = new Screen(new Vector2(-8, -8), new Vector2(8, 8), width, height, 1);
        PerspectiveProjector projector = new PerspectiveProjector(5);

        JFrame frame = new JFrame();
        frame.setSize(width, height);
        frame.getContentPane().setLayout(new FlowLayout());
        JLabel label = new JLabel(new ImageIcon(screen.image()));
        frame.getContentPane().add(label);
        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        long t0 = System.nanoTime();
        long t1 = System.nanoTime();

        boolean move = true;
        double speedx = 0;
        double amplitudex = 0;
        double speedy = 2;
        double amplitudey = 1.5;
        double speedz = 0;
        double amplitudez = 0;
        double rotation_speedx = 0.5;
        double rotation_speedy = 1.5;
        //double speed_translation = 2;

        //int nCube = 3;
        //Cube[] cubes = new Cube[nCube];
        //for(int i = 0; i < nCube; ++i) {
        //    cubes[i] = Cube.randomCube(random, -2, 2, -2, 2, 2, 4, 0.8, 2.0);
        //}
        //int iter = 0;
        //int replacement_speed = 1;
        //int replacements = 0;
        //int replacement_wait = 15;

        //int cube_square_size = 30;
        //int nCubes = cube_square_size * cube_square_size;
        //double minx_cubes = -1;
        //double maxx_cubes = 1;
        //double minz_cubes = -1;
        //double maxz_cubes = 5;
        //Object3D[] objects = new Object3D[nCubes];
        //for(int i = 0; i < nCubes; ++i) {
        //    objects[i] = new Cube(
        //        new Vector3((double)(i % cube_square_size) * (maxx_cubes - minx_cubes) / (cube_square_size - 1) + minx_cubes,
        //            new Random().nextDouble() * 4 - 2, 
        //                    (double) (i / cube_square_size) * (maxz_cubes - minz_cubes) / (cube_square_size - 1) + minz_cubes
        //                ), 
        //             (maxx_cubes - minx_cubes) / cube_square_size, 0xffff00, i * 8);
        //}
        //Cube cube1 = new Cube(new Vector3(-2, -1.5, 4), 2, 0x77ff00ff, 0);
        //Vector3 cube_center1 = new Vector3(-1, -0.5, 5);
        //Cube cube2 = new Cube(new Vector3(1, -1.5, 4), 2, 0x77ff00ff, 8);
        //Vector3 cube_center2 = new Vector3(2, -0.5, 5);

        Object3D[] objects = { Object3D.readObj(new File("models/teapot.obj")) };

        Vector3[] base_points = get_points(objects);
        Triangle triangles[] = get_triangles(objects);
        for (Triangle t : triangles) {
            t.compute_inner_triangle_normal(base_points);
        }

        Vector3[] light_dirs = {
            new Vector3(1, -1, 1).normalized(),
            new Vector3(-1, -1, 1.2).normalized(),
            //new Vector3(0, -1, 0).normalized()
        };

        while(true) {

            //time stuff
            try {
                Thread.sleep(1);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            long current_time = System.nanoTime();
            System.out.printf("Frame time : %f ms\n", (current_time - t1) * 1e-6);
            t1 = current_time;
            double elapsed = (t1 - t0) * 1e-9;


            //Vector3 offset = move ? new Vector3(Math.sin((t1 - t0)* 1e-9 * speedx) * amplitudex, 
            //    Math.sin((t1 - t0)* 1e-9 * speedy) * amplitudey, 0) : new Vector3(0,0,0);
            //Vector3 offset = new Vector3(0, 0, - elapsed * speed_translation);

            //Cube cube = new Cube(new Vector3(-1, -1, 2).add(offset), 1);
            //Object3D[] objects = new Object3D[nCube];
            //for(int i = 0; i < nCube; ++i) {
            //    objects[i] = cubes[i].transform(offset, 1);
            //}
            //if(replacements % replacement_wait == 0) {
            //    for(int i = 0; i < replacement_speed; ++i) {
            //        cubes[iter % nCube] = Cube.randomCube(random, -2, 2, -2, 2, elapsed * speed_translation + 2, elapsed * speed_translation + 20, 0.8, 2);
            //        iter++;
            //    }
            //}
            //replacements++;

            //Project them onto 2D plane with a perspective projector
            //Vector2[][] object_projected_points = new Vector2[objects.length][];
            //for(int i = 0; i < objects.length; ++i) {
            //    Vector2[] projected_points = new Vector2[objects[i].points().length];
            //    for(int j = 0; j < objects[i].points().length; ++j) {
            //        projected_points[j] = objects[i].points()[j].z() < 0.5 ? null : projector.project(objects[i].points()[j]);
            //    }

            //    object_projected_points[i] = projected_points;
            //}

            long tbt = System.nanoTime();

            Vector3 offset = move ? new Vector3(Math.sin((t1 - t0)* 1e-9 * speedx) * amplitudex, 
                Math.sin((t1 - t0)* 1e-9 * speedy) * amplitudey, Math.sin((t1 - t0)* 1e-9 * speedz) * amplitudez + 5) : new Vector3(2,2, 2);
            
            int n_threads = 8;
            Thread[] threads = new Thread[n_threads];
            double ratio = (double) base_points.length / n_threads;
            Vector3[] points = new Vector3[base_points.length];
            for(int k = 0; k < n_threads; ++k) {
                final int start = (int) (ratio * k);
                final int end = (int) (ratio * (k + 1));
                final long ft1 = t1;
                threads[k] = new Thread(() -> {
                    for(int i = start; i < end; ++i) {
                        points[i] = base_points[i].plus(offset);
                        points[i] = points[i].rotate_y_axis(offset, (ft1 - t0) * 1e-9 * rotation_speedy);
                        points[i] = points[i].rotate_x_axis(offset, (ft1 - t0) * 1e-9 * rotation_speedx);
                    }
                });
                threads[k].start();
            }
            for(Thread t : threads) {
                try {
                    t.join();
                } catch (InterruptedException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }

            long tat = System.nanoTime();
            double elapsed_transform = (tat - tbt) * 1e-6;
            System.out.println("Time to transform : " + elapsed_transform + "ms");

            long tbd = System.nanoTime();
            //Draw the plane by defining a screen on it
            screen.drawTriangles(triangles, points, projector, light_dirs);
            long tad = System.nanoTime();
            double elapsed_draw = (tad - tbd) * 1e-6;
            System.out.println("Time to draw : " + elapsed_draw + "ms");

            long tbr = System.nanoTime();
            label.setIcon(new ImageIcon(screen.image()));

            screen.resetImage();
            long tar = System.nanoTime();
            double elapsed_reset = (tar - tbr) * 1e-6;
            System.out.println("Time to reset : " + elapsed_reset + "ms");

        }

        //screen.save("png", new File("result.png"));

    }

    private static Vector3[] get_points(Object3D[] objects) {

        int length = 0;
        for(Object3D o : objects) length += o.points().length;

        Vector3[] result = new Vector3[length];

        int i = 0;
        for(Object3D o : objects) {
            for (Vector3 v : o.points()) {
                result[i] = v;
                i += 1;
            }
        }

        return result;

    }

    private static Triangle[] get_triangles(Object3D[] objects) {

        int length = 0;
        for(Object3D o : objects) length += o.triangles().length;

        Triangle[] result = new Triangle[length];

        int i = 0;
        for(Object3D o : objects) {
            for (Triangle v : o.triangles()) {
                result[i] = v;
                i += 1;
            }
        }

        return result;

    }

}

interface Object3D {
    
    Vector3[] points();

    Triangle[] triangles();

    int color();

    public static Object3D readObj(File file) {

        List<String> lines = new ArrayList<>();

        try (BufferedReader bf = new BufferedReader(new FileReader(file))) {
            String l;
            while((l = bf.readLine()) != null) {
                lines.add(l);
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        int nVertices = 0;
        int nTriangles = 0;
        for(String s : lines) {
            if(s.startsWith("v ")) nVertices++;
            if(s.startsWith("f ")) nTriangles++;
        }

        Vector3[] vertices = new Vector3[nVertices];
        Triangle[] triangles = new Triangle[nTriangles]; 

        int iV = 0;
        int iT = 0;
        for(String s : lines) {

            String[] split = s.split(" ");

            if(split[0].equals("v")) {
                vertices[iV] = new Vector3(Double.parseDouble(split[1]), Double.parseDouble(split[2]), Double.parseDouble(split[3]));
                iV += 1;
            } else if(split[0].equals("f")) {
                triangles[iT] = new Triangle(Integer.parseInt(split[1]) - 1, Integer.parseInt(split[2]) - 1, Integer.parseInt(split[3]) - 1);
                iT += 1;
            }
        }

        return new GeneralObject3D(vertices, triangles);

    } 

}

class GeneralObject3D implements Object3D {

    Vector3[] points;
    Triangle[] triangles;

    public GeneralObject3D(Vector3[] points, Triangle[] triangles) {
        this.points = points;
        this.triangles = triangles;
    }

    public Vector3[] points() {
        return points;
    }

    public Triangle[] triangles() {
        return triangles;
    }

    public int color() {
        return 0;
    }
 
}

class Cube implements Object3D {

    private final Vector3 bottomLeft;
    private final double size;

    private final int color;

    private final Vector3[] points;

    private final Triangle[] triangles;
    private final int offset;

    public static int[][] cube_triangle_data = {
        {0, 1, 3},
        {0, 3, 2},
        {4, 7, 5},
        {4, 6, 7},
        {0, 5, 1},
        {0, 4, 5},
        {0, 6, 4},
        {0, 2, 6},
        {2, 3, 7},
        {2, 7, 6},
        {1, 5, 7},
        {1, 7, 3}
    };

    public Cube(Vector3 bottomLeft, double size, int color, int offset) {
        this.bottomLeft = bottomLeft;
        this.size = size;

        this.color = color;

        this.points = new Vector3[8];
        
        for(int i = 0; i < 2; ++i) {
            for(int j = 0; j < 2; ++j) {
                for(int k = 0; k < 2; ++k) {
                    points[i + 2*j + 4*k] = bottomLeft.plus(new Vector3(i, j, k).scale(size));
                }
            }
        }

        this.offset = offset;
        this.triangles = new Triangle[12];
        for(int i = 0; i < triangles.length; ++i) {
            triangles[i] = new Triangle(cube_triangle_data[i][0] + offset, cube_triangle_data[i][1] + offset, cube_triangle_data[i][2] + offset);
        }

    }

    public static Cube randomCube(Random random, double x1, double x2, double y1, double y2, double z1, double z2, double s1, double s2, int offset) {

        double x = random.nextDouble() * (x2 - x1) + x1;
        double y = random.nextDouble() * (y2 - y1) + y1;
        double z = random.nextDouble() * (z2 - z1) + z1;
        double s = random.nextDouble() * (s2 - s1) + s1;

        return new Cube(new Vector3(x, y, z), s, 0x44cc99, offset);
    } 

    public Cube transform(Vector3 offset, double scale, int i_offset) {
        return new Cube(this.bottomLeft.plus(offset), this.size * scale, this.color, i_offset);
    }

    public Vector3[] points() {
        return points;
    }

    public Triangle[] triangles() {
        return triangles;
    }

    public int color() {
        return color;
    }
}

class Vector2{

    final private double x;
    final private double y;

    public Vector2(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double norm() {
        return Math.sqrt(x*x + y*y);
    }

    public Vector2 plus(Vector2 that) {
        return new Vector2(this.x + that.x, this.y + that.y);
    }

    public Vector2 minus(Vector2 that) {
        return new Vector2(this.x - that.x, this.y - that.y);
    }

    public Vector2 scale(double scale) {
        return new Vector2(this.x * scale, this.y * scale);
    }

    public double x() {
        return x;
    }

    public double y() {
        return y;
    }

    public Vector2 normalized() {
        return this.scale(1 / this.norm());
    }

    public double dot(Vector2 that) {
        return this.x * that.x + this.y * that.y;
    }

    public double cross(Vector2 that) {
        return this.x * that.y - this.y * that.x;
    }

}

class Vector3{

    final private double x;
    final private double y;
    final private double z;

    public Vector3(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double norm() {
        return Math.sqrt(x*x + y*y + z*z);
    }

    public Vector3 plus(Vector3 that) {
        return new Vector3(this.x + that.x, this.y + that.y, this.z + that.z);
    }

    public Vector3 minus(Vector3 that) {
        return new Vector3(this.x - that.x, this.y - that.y, this.z - that.z);
    }

    public Vector3 scale(double scale) {
        return new Vector3(this.x * scale, this.y * scale, this.z * scale);
    }

    public double x() {
        return x;
    }

    public double y() {
        return y;
    }

    public double z() {
        return z;
    }

    public Vector3 normalized() {
        return this.scale(1 / this.norm());
    }

    public double dot(Vector3 that) {
        return this.x * that.x + this.y * that.y + this.z * that.z;
    }

    public Vector3 cross(Vector3 that) {
        return new Vector3(
            this.y * that.z - this.z * that.y,
            this.z * that.x - this.x * that.z,
            this.x * that.y - this.y * that.x
        );
    }

    public Vector3 rotate_y_axis(Vector3 center, double angle) {

        Vector3 centered = this.minus(center);

        double cos = Math.cos(angle);
        double sin = Math.sin(angle);
        Vector3 rotated = new Vector3(
            cos * centered.x + sin * centered.z,
            centered.y,
            -sin * centered.x + cos * centered.z
        );

        return rotated.plus(center);
    }

    public Vector3 rotate_x_axis(Vector3 center, double angle) {

        Vector3 centered = this.minus(center);

        double cos = Math.cos(angle);
        double sin = Math.sin(angle);
        Vector3 rotated = new Vector3(
            centered.x,
            cos * centered.y - sin * centered.z,
            sin * centered.y + cos * centered.z
        );

        return rotated.plus(center);
    }

    @Override
    public String toString() {
        return String.format("[%f, %f, %f]", x, y, z);
    }

}

class PerspectiveProjector {

    private final double f;
    
    public PerspectiveProjector(double f) {
        this.f = f;
    }

    public Vector2 project(Vector3 v) {

        double mu = f / v.z();

        return new Vector2(mu * v.x(), mu * v.y());
    }
}

class Circle {

    private final Vector2 pos;
    private final double r;

    public Circle(Vector2 pos, double r) {
        this.pos = pos;
        this.r = r;
    }

    public boolean contains(Vector2 p) {
        return pos.plus(p.scale(-1)).norm() < r;
    }

    public Vector2 center() {
        return pos;
    }

}

class Sphere {

    private final Vector3 pos;
    private final double r;

    public Sphere(Vector3 pos, double r) {
        this.pos = pos;
        this.r = r;
    }

    public boolean contains(Vector3 p) {
        return pos.plus(p.scale(-1)).norm() < r;
    }

    public Circle toCircle() {
        return new Circle(new Vector2(pos.x(), pos.y()), r);
    }

    public Vector3 center() {
        return pos;
    }

}

class Triangle {

    public final int i1, i2, i3;

    public final double r = 0.871; 
    public final double g = 0.667; 
    public final double b = 0.647; 


    //cached normal vectors to the sides of the triangle
    public boolean computed = false;
    private Vector2 n1, n2, n3;

    //cached 3D inner normal of the triangle
    private Vector3 n;

    public Triangle(int i1, int i2, int i3) {

        this.i1 = i1;
        this.i2 = i2; 
        this.i3 = i3;

    }

    public boolean contains(Vector2 point, Vector2[] points) {

        Vector2 p1 = points[i1], p2 = points[i2], p3 = points[i3];
        if(p1 == null || p2 == null || p3 == null) return false;

        //for each side, find the line equation that is positive if the pixel is on the correct side, and return the conjuction of the 3 signs.
        //equation is of the form ax + by + c = 0, or n * (p - p1) = 0 for a side containing p1 for example.
        if(!computed) {
            n1 = inner_side_normal(p1, p2, p3);
            n2 = inner_side_normal(p2, p3, p1);
            n3 = inner_side_normal(p3, p1, p2);

            computed = true;
        }

        boolean b1 = n1.dot(point.minus(p1)) >= 0;
        boolean b2 = n2.dot(point.minus(p2)) >= 0;
        boolean b3 = n3.dot(point.minus(p3)) >= 0;

        return b1 && b2 && b3;

    }

    public Vector3 p1(Vector3[] points) {return points[i1];}
    public Vector3 p2(Vector3[] points) {return points[i2];}
    public Vector3 p3(Vector3[] points) {return points[i3];}

    public Vector2 p1(Vector2[] points) {return points[i1];}
    public Vector2 p2(Vector2[] points) {return points[i2];}
    public Vector2 p3(Vector2[] points) {return points[i3];}

    public void compute_inner_triangle_normal(Vector3[] points) {

        //Goal is to find a normal vector pointing to the interior of the object the triangle is in

        Vector3 p1 = points[i1], p2 = points[i2], p3 = points[i3];

        Vector3 d1 = p2.minus(p1), d2 = p3.minus(p1);

        n = d1.cross(d2).normalized().scale(-1);
    }

    public Vector3 normal() {
        return n;
    }

    private Vector2 inner_side_normal(Vector2 p1, Vector2 p2, Vector2 p3) {

        //here we look for a normal vector to the line defined by p1 and p2, extending in the semi-plan where p3 is.
        //let's consider d = p2 - p1 = (p2.x - p1.x, p2.y - p1.y). What we want is that n * d = 0 basically.
        //so n.x * (p2.x - p1.x) + n.y * (p2.y - p1.y) = 0. If I isolate n.x : 
        //n.x = ( -n.y * (p2.y - p1.y) ) / ( p2.x - p1.x )
        //we kind of have a system, variables are n.x and n.y. another constraint : n * (p3 - p1) > 0

        if(p1.x() == p2.x()) {

            //they have the same x, normal must be on x only
            return new Vector2(Math.signum(p3.x() - p1.x()), 0);

        } else {

            double ny = 1;
            double nx = ( -ny * (p2.y() - p1.y()) ) / ( p2.x() - p1.x() );

            Vector2 n = new Vector2(nx, ny);

            return n.scale(Math.signum(n.dot(p3.minus(p1))));

        }

    }

    @Override
    public String toString() {
        return String.format("Triangle : \n   Ids : (%d, %d, %d)\n   Normal : " + n.toString(), i1, i2, i3);
    }

}

class Line {

    private final int i1, i2;

    public Line(int i1, int i2) {
        this.i1 = i1;
        this.i2 = i2;
    }

    public Vector2 start(Vector2[] points) {
        return points[i1];
    }

    public Vector2 end(Vector2[] points) {
        return points[i2];
    }

    public Vector2 dir(Vector2[] points) {

        Vector2 diff = end(points).plus(start(points).scale(-1));

        return diff.normalized();

    }

    public double length(Vector2 points[]) {
        return end(points).plus(start(points).scale(-1)).norm();
    }

}

class Screen {

    private final Vector2 bottomLeft, topRight;
    private final double rlx, rly;
    private final double rli, rlj;

    private final int width, height;
    private final int upscale;
    private BufferedImage image;

    public Comparator<Triangle> zTriComparator(Vector3[] points) {
        return (t1, t2) -> {
            double z1 = Math.min(t1.p1(points).z(), Math.min(t1.p2(points).z(), t1.p3(points).z())); 
            double z2 = Math.min(t2.p1(points).z(), Math.min(t2.p2(points).z(), t2.p3(points).z())); 

            return Double.compare(z2, z1);
        };
    }

    public Screen(Vector2 bottomLeft, Vector2 topRight, int width, int height, int upscale) {
        
        this.bottomLeft = bottomLeft;
        this.topRight = topRight;
        this.width =  width;
        this.height = height;
        this.upscale = upscale;

        this.rlx = (topRight.x() - bottomLeft.x()) / width;
        this.rly = (topRight.y() - bottomLeft.y()) / height;
        this.rli = 1 / this.rlx;
        this.rlj = 1 / this.rly;

        this.image = new BufferedImage(width * upscale, height * upscale, BufferedImage.TYPE_INT_RGB);
    }

    private Vector2 point_from_pixel(int i, int j) {
        return new Vector2( rlx * i + x1(), rly * j + y1());
    }

    private int i_from_x(double x) {
        return (int) (rli * (x - x1()));
    }

    private int j_from_y(double y) {
        return (int) (rlj * (y - y1()));
    }

    public void drawPoint(Vector2 point, int color) {

        int i = i_from_x(point.x());
        int j = j_from_y(point.y());

        drawPixel(i, j, color);        

    }

    public void drawPixel(int i, int j, int color) {
        if(i >= 0 && i < width && j >= 0 && j < height) {
            //safe to draw
            for(int k = 0; k < upscale; ++k) {
                for(int l = 0; l < upscale; ++l) {
                    image.setRGB(i * upscale + k, (height - j - 1) * upscale + l, color);
                }  
            }
        }
    }

    public void drawLine(Line line, Vector2[] points) {
        //System.out.println("heere");
        //we need to know what 1 distance in pixel place is in plane space
        if(line.end(points) == null || line.start(points) == null) return;

        double d = (x2() - x1()) / width;
        Vector2 dv = line.dir(points).scale(d);
        
        Vector2 currentPos = line.start(points);
        int steps = (int) (line.length(points) / d);
        //System.out.println("d is " + d + ", steps is " + steps);
        for(int i = 0; i < steps; ++i) {
            //System.out.println("there, " + i);
            drawPoint(currentPos, 0xaaffaa);
            currentPos = currentPos.plus(dv);
        }
    }

    private double illumination(Triangle t, Vector3 light_dir) {

        double k_d = 0.4;
        double k_s = 1.0;
        double ambiant = 0.1;
        double cst = 0.8;

        double diffuse = Math.max(0, light_dir.dot(t.normal()));
        diffuse *= k_d;

        Vector3 v = new Vector3(0, 0, 1);
        Vector3 h = light_dir.plus(v).normalized();
        double q = 70;
        double specular = k_s * Math.pow(t.normal().dot(h), q);

        return ambiant + (diffuse + specular) * cst;
    }

    private void rasterize(Triangle[] triangles, int start, int end, Vector3[] points, Vector2[] projected_points, VertexShadingState[] vSS) {
        for(int k = start; k < end; ++k) {
            Triangle t = triangles[k];
            if(t.p1(projected_points) == null || t.p2(projected_points) == null || t.p3(projected_points) == null) continue;

            if(t.normal().dot(t.p1(points)) <= 0) continue;

            //find min and max indices for i and j

            double minx = Math.min(Math.min(t.p1(projected_points).x(), t.p2(projected_points).x()), t.p3(projected_points).x());
            double miny = Math.min(Math.min(t.p1(projected_points).y(), t.p2(projected_points).y()), t.p3(projected_points).y());
            double maxx = Math.max(Math.max(t.p1(projected_points).x(), t.p2(projected_points).x()), t.p3(projected_points).x());
            double maxy = Math.max(Math.max(t.p1(projected_points).y(), t.p2(projected_points).y()), t.p3(projected_points).y());

            int mini = Math.max(1, i_from_x(minx)); //clamp to w and h for performance
            int minj = Math.max(1, j_from_y(miny));
            int maxi = Math.min(width-1, i_from_x(maxx));
            int maxj = Math.min(height-1, j_from_y(maxy));

            Vector3 average_pos = t.p1(points).plus(t.p2(points)).plus(t.p3(points));

            for(int i = mini-1; i < maxi+1; ++i) {
                //boolean previous_contained = false;
                for(int j = minj-1; j < maxj+1; ++j) {
                    Vector2 point = point_from_pixel(i, j);

                    if(t.contains(point, projected_points)) {
                        //if(previous_contained) {
                            //Shading stuff

                            double d1, d2, d3; //the inverted distances to the 3 points.
                            d1 = 1 / point.minus(t.p1(projected_points)).norm();
                            d2 = 1 / point.minus(t.p2(projected_points)).norm();
                            d3 = 1 / point.minus(t.p3(projected_points)).norm();

                            double sum = d1 + d2 + d3;
                            double t1 = d1 / sum;
                            double t2 = d2 / sum;
                            double t3 = d3 / sum;

                            double interp = t1 * vSS[t.i1].intensity + t2 * vSS[t.i2].intensity + t3 * vSS[t.i3].intensity;

                            double r = Math.sin(average_pos.x() / 10.0) * 0.5 + 0.5;
                            double g = Math.sin(average_pos.y() / 10.0) * 0.5 + 0.5;
                            double b = Math.sin(average_pos.z() / 10.0) * 0.5 + 0.5;

                            r = Math.min(r * interp, 1);
                            g = Math.min(g * interp, 1);
                            b = Math.min(b * interp, 1);

                            int ir = (int) (r * 0xff);
                            int ig = (int) (g * 0xff);
                            int ib = (int) (b * 0xff);

                            int rgb = (ir << 16) | (ig << 8) | (ib);

                            drawPixel(i, j, rgb);
                         /* } else {
                            previous_contained = true;
                            drawPixel(i, j, 0xffffff);
                        }
                     } else {
                        if(previous_contained) {
                            drawPixel(i, j, 0xffffff);
                        }
                        previous_contained = false;*/
                    }
                }
            }
        }
    }

    class VertexShadingState {
        double intensity = 0;
        int triangle_count = 0;
    }

    public void drawTriangles(Triangle[] triangles, Vector3[] points, PerspectiveProjector projector, Vector3[] light_dirs) {

        long tbp = System.nanoTime();
        Vector2[] projected_points = new Vector2[points.length];
        for(int i = 0; i < points.length; ++i) {
            projected_points[i] = points[i].z() < 0.2  ? null : projector.project(points[i]);
        }

        VertexShadingState[] vSS = new VertexShadingState[points.length];
        for(Triangle t : triangles) {
            t.computed = false;
            t.compute_inner_triangle_normal(points);

            //illuminate vertices
            double ill = illumination(t, light_dirs[0]);
            for(int i : List.of(t.i1, t.i2, t.i3)) {
                if(vSS[i] == null) {
                    vSS[i] = new VertexShadingState();
                    vSS[i].intensity = ill;
                    vSS[i].triangle_count = 1;           
                } else {
                    vSS[i].intensity = vSS[i].intensity * vSS[i].triangle_count / (vSS[i].triangle_count + 1) + ill / (vSS[i].triangle_count + 1);
                    vSS[i].triangle_count++;
                }
            }
        }

        Arrays.sort(triangles, zTriComparator(points));
        long  tap = System.nanoTime();
        double elasped_projection = (tap - tbp) * 1e-6;
        System.out.println("Time to project : " + elasped_projection + "ms");

        int n_threads = 1;
        Thread[] threads = new Thread[n_threads];
        long t0 = System.nanoTime();
        double ratio = (double) triangles.length / n_threads; //the amount of triangles to handle per thread
        for(int i = 0; i < n_threads; ++i) {
            final int start = (int) (ratio * i);
            final int end = (int) (ratio * (i+1));
            threads[i] = new Thread(() -> {
                rasterize(triangles, start, end, points, projected_points, vSS);
            });
            threads[i].start();
        }

        for(Thread t : threads) {
            try {
                t.join();
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        double dt = (System.nanoTime() - t0) * 1e-6;
        System.out.printf("Time to raster : %fms \n", dt);

    }


    public void save(String format, File file) {

        try {
            ImageIO.write(image, format, file);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

    public BufferedImage image() {
        return image;
    }

    public void resetImage() {
        this.image = new BufferedImage(width*upscale, height*upscale, image.getType());
        return;
    }

    public Vector2 bottomLeft() {
        return bottomLeft;
    }

    public Vector2 topRight() {
        return topRight;
    }

    public double x1() {
        return bottomLeft.x();
    }

    public double y1() {
        return bottomLeft.y();
    }

    public double x2() {
        return topRight.x();
    }

    public double y2() {
        return topRight.y();
    }

}
