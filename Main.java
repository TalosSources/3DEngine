import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Random;

import javax.imageio.ImageIO;
import java.io.File;

import javax.swing.*;
import java.awt.*;

public class Main{

    private static final Random random = new Random();

    public static void main(String[] args) {

        Screen screen = new Screen(new Vector2(-8, -8 * 108.0 / 192.0), new Vector2(8, 8 * 108.0 / 192.0), 576, 324, 3);
        PerspectiveProjector projector = new PerspectiveProjector(5);

        JFrame frame = new JFrame();
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
        double speedy = 3;
        double amplitudey = 0.7;
        double rotation_speed = 2;
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

        Cube cube = new Cube(new Vector3(-1, -1.5, 4), 2, 0x77ff00ff);
        Vector3 cube_center = new Vector3(0, -0.5, 5);

        Vector3[] base_points = cube.points();
        Triangle triangles[] = cube.triangles();
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
            t1 = System.nanoTime();
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

            Vector3 offset = move ? new Vector3(Math.sin((t1 - t0)* 1e-9 * speedx) * amplitudex, 
                Math.sin((t1 - t0)* 1e-9 * speedy) * amplitudey, 0) : new Vector3(2,2, 2);
            

            Vector3[] points = new Vector3[base_points.length];
            for(int i = 0; i < points.length; ++i) {
                points[i] = base_points[i].plus(offset);
                points[i] = points[i].rotate_y_axis(cube_center.plus(offset), (t1 - t0) * 1e-9 * rotation_speed);
            }

            screen.resetImage();

            //Draw the plane by defining a screen on it
            screen.drawTriangles(triangles, points, projector, light_dirs);

            label.setIcon(new ImageIcon(screen.image()));

        }

        //screen.save("png", new File("result.png"));

    }
}

interface Object3D {
    
    Vector3[] points();

    Triangle[] triangles();

    int color();

}

class Cube implements Object3D {

    private final Vector3 bottomLeft;
    private final double size;

    private final int color;

    private final Vector3[] points;
    //private final Line[] lines = {
    //    new Line(0, 1),
    //    new Line(1, 3),
    //    new Line(3, 2),
    //    new Line(2, 0),
    //    new Line(4, 5),
    //    new Line(5, 7),
    //    new Line(7, 6),
    //    new Line(6, 4),
    //    new Line(0, 4),
    //    new Line(1, 5),
    //    new Line(2, 6),
    //    new Line(3, 7)
    //};

    private final Triangle[] triangles = {
        new Triangle(0, 1, 3, 4), //front
        new Triangle(0, 2, 3, 4),
        new Triangle(4, 5, 7, 0), //back
        new Triangle(4, 6, 7, 0),
        new Triangle(0, 1, 5, 7), //bottom
        new Triangle(0, 4, 5, 7),
        new Triangle(0, 4, 6, 7), //left
        new Triangle(0, 2, 6, 7),
        new Triangle(2, 3, 7, 0), //top
        new Triangle(2, 6, 7, 0),
        new Triangle(1, 5, 7, 0), //right
        new Triangle(1, 3, 7, 0)
    };

    public Cube(Vector3 bottomLeft, double size, int color) {
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

    }

    public static Cube randomCube(Random random, double x1, double x2, double y1, double y2, double z1, double z2, double s1, double s2) {

        double x = random.nextDouble() * (x2 - x1) + x1;
        double y = random.nextDouble() * (y2 - y1) + y1;
        double z = random.nextDouble() * (z2 - z1) + z1;
        double s = random.nextDouble() * (s2 - s1) + s1;

        return new Cube(new Vector3(x, y, z), s, random.nextInt());
    } 

    public Cube transform(Vector3 offset, double scale) {
        return new Cube(this.bottomLeft.plus(offset), this.size * scale, this.color);
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
            y,
            -sin * centered.x + cos * centered.z
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

    private final int i1, i2, i3;
    private final int ref_id;

    public final double r = new Random().nextDouble(); 
    public final double g = new Random().nextDouble(); 
    public final double b = new Random().nextDouble(); 


    //cached normal vectors to the sides of the triangle
    public boolean computed = false;
    private Vector2 n1, n2, n3;

    //cached 3D inner normal of the triangle
    private Vector3 n;

    public Triangle(int i1, int i2, int i3, int ref_id) {

        this.i1 = i1;
        this.i2 = i2; 
        this.i3 = i3;

        this.ref_id = ref_id;

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

    public void compute_inner_triangle_normal(Vector3[] points) {

        //Goal is to find a normal vector pointing to the interior of the object the triangle is in

        Vector3 p1 = points[i1], p2 = points[i2], p3 = points[i3], ref = points[ref_id];

        Vector3 d1 = p2.minus(p1), d2 = p3.minus(p1);

        n = d1.cross(d2);

        n = n.scale(Math.signum(n.dot(ref.minus(p1)))).normalized();

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
    private boolean print = true;

    private final Vector2 bottomLeft, topRight;
    private final int width, height;
    private final int upscale;
    private final BufferedImage image;

    public Screen(Vector2 bottomLeft, Vector2 topRight, int width, int height, int upscale) {
        
        this.bottomLeft = bottomLeft;
        this.topRight = topRight;
        this.width =  width;
        this.height = height;
        this.upscale = upscale;

        this.image = new BufferedImage(width * upscale, height * upscale, BufferedImage.TYPE_INT_RGB);
    }

    public void drawPoint(Vector2 point, int color) {

        int i = (int) (width * (point.x() - x1()) / (x2() - x1()));
        int j = (int) (height * (topRight.y() - point.y() - y1()) / (y2() - y1()));

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

    public void drawTriangles(Triangle[] triangles, Vector3[] points, PerspectiveProjector projector, Vector3[] light_dirs) {

        Vector2[] projected_points = new Vector2[points.length];
        for(int i = 0; i < points.length; ++i) {
            projected_points[i] = points[i].z() < 0.5  ? null : projector.project(points[i]);
        }
        for(Triangle t : triangles) {
            t.computed = false;
            t.compute_inner_triangle_normal(points);
        }

        for(int i = 0; i < width; ++i) {
            for(int j = 0; j < height; ++j) {

                Vector2 point = new Vector2((x2() - x1()) * i / width + x1(), (y2() - y1()) * j / height + y1());

                for(Triangle t : triangles) {
                    if(print)
                        System.out.println("T is " + t + ", p1 is " + t.p1(points));
                    if(t.normal().dot(t.p1(points)) <= 0) continue;
                    if(t.contains(point, projected_points)) {
                        //Shading stuff

                        double r = t.r;
                        double g = t.g;
                        double b = t.b;

                        double intensity = 0;
                        for (Vector3 light_dir : light_dirs) {
                            intensity += Math.max(0, light_dir.dot(t.normal()));
                        } 
                        //intensity = 1.0;
                        r = Math.min(r * intensity, 1);
                        g = Math.min(g * intensity, 1);
                        b = Math.min(b * intensity, 1);

                        int ir = (int) (r * 0xff);
                        int ig = (int) (g * 0xff);
                        int ib = (int) (b * 0xff);

                        int rgb = (ir << 16) | (ig << 8) | (ib);

                        drawPixel(i, j, rgb);

                        break;

                    }
                }
                print = false;
            }
        }

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
        for(int i = 0; i < width*upscale; ++i) {
            for(int j = 0; j < height*upscale; ++j) {
                image.setRGB(i, j, 0xff000000);
            } 
        }
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