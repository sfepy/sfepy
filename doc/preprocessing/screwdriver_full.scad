radius = 0.01;
height = 0.1;
$fn = 50;

module tip() {
  rotate([0, -10, 0])
    translate([0, -radius, -3*radius])
      cube([radius, 2*radius, 3*radius], center=false);
}

difference() {
  difference() {
    difference() {
      union() {
        cylinder(center=false, h=height, r=radius);
        sphere(radius);
      };
      translate([0, 0, 0.9*height])
        rotate_extrude()
          polygon([[0.8*radius, 0], [1.8*radius, -0.577*radius], [1.8*radius, 0.577*radius]]);
    }
    cylinder(center=false, h=height, r=0.3*radius);
  }
  for (i = [1:6]) {
    rotate([0, 0, 360/6*i])
      translate([-1.1*radius, 0.0, -0.2*height])
        cylinder(center=false, h=1.1*height, r=0.2*radius);
  }
}

union() {
  difference() {
    translate([0, 0, height])
      cylinder(center=false, h=height, r=0.3*radius);
    translate([0, 0, 1.71*height + 3*radius])
      union() {
        tip();
        mirror ([1, 0, 0]) tip();
      }
  }
  cylinder(center=false, h=height, r=0.3*radius);
}
