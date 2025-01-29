class Point {
	constructor(id = null, x = null, y = null, attributes = null, boundaryMarker = null) {
		this.Id = id;
		this.X = x;
		this.Y = y;
		this.Attributes = attributes;
		this.BoundaryMarker = boundaryMarker;
	}

	equals(other) {
		return (this.X === other.X && this.Y === other.Y);
	}
}

class Vertex extends Point {
	constructor(id = null, x = null, y = null, attributes = null, boundaryMarker = null) {
		super(id, x, y, attributes, boundaryMarker);
		this.Mark = 0;
		this.Hash = 0;
		this.Tri = new Otri();
		this.Type = VertexType.InputVertex;
	}
}

// Pointer to a region in the mesh geometry. A region is a well-defined
// subset of the geomerty (enclosed by subsegments).
class RegionPointer {
	constructor(x, y, id) {
		this.Id = id;
		this.Point = new Point(-1, x, y);
	}
}

class Hole {
	constructor(id, x, y) {
		this.Id = id;
		this.X = x;
		this.Y = y;
	}
}
class Behavior {
	constructor() {
		this.Poly = false;
		this.VarArea = false;
		this.UseRegions = false;
		this.SteinerPoints = -1;
		this.UseSegments = true;
		this.Convex = false;
	}
}

class Edge {
	constructor(p0, p1, boundary = 0) {
		this.P0 = p0;
		this.P1 = p1;
		this.Boundary = boundary;
	}
}
class Sampler {
	constructor() {
		// Number of random samples for point location (at least 1).
		this.Samples = 1;
		// Number of triangles in mesh.
		this.TriangleCount = 0;
		// Empirically chosen factor.
		this.Samplefactor = 11;
		// Keys of the triangle dictionary.
		this.Keys = [];
	}

	// Reset the sampler.
	reset() {
		this.Samples = 1;
		this.TriangleCount = 0;
	}

	//Update sampling parameters if mesh changed.
	update(mesh) {
		this.updateMesh(mesh, false);
	}

	//Update sampling parameters if mesh changed.
	updateMesh(mesh, forceUpdate) {
		let count = mesh.Triangles.size;
		if (this.TriangleCount !== count || forceUpdate) {
			this.TriangleCount = count;
			// The number of random samples taken is proportional to the cube root of
			// the number of triangles in the mesh.
			while (this.Samplefactor * this.Samples * this.Samples * this.Samples < count) {
				this.Samples++;
			}
			this.Keys = Array.from(mesh.Triangles.keys());
		}
	}

	getSamples(mesh) {
		let randSamples = [];
		let range = this.TriangleCount / this.Samples;
		let key = 0;
		for (var i = 0; i < this.Samples; i++) {
			let min = Math.ceil(i * range);
			let max = Math.floor((i + 1) * range - 1);
			let key = Math.floor(Math.random() * (max - min)) + min;
			if (!mesh.Triangles.has(this.Keys[key])) {
				this.updateMesh(mesh, true);
				i--;
			} else {
				randSamples.push(this.Keys[key]);
			}
		}
		return randSamples;
	}

}

class InputGeometry {
	constructor(numOfVertices,
		dimension,
		verticesNumOfAttributes,
		verticesNumOfBoundaryMakers,
		vertices,
		boundingBox,
		numOfSegments,
		segmentsNumOfBoundaryMakers,
		numberOfHoles,
		numberOfRegionalAttrs,
		segments,
		holes) {

		this.NumOfVertices = numOfVertices;
		this.Dimension = dimension;
		this.VerticesNumOfAttributes = verticesNumOfAttributes;
		this.VerticesNumOfBoundaryMakers = verticesNumOfBoundaryMakers;
		this.Vertices = vertices;
		this.BoundingBox = boundingBox;
		this.NumOfSegments = numOfSegments;
		this.SegmentsNumOfBoundaryMakers = segmentsNumOfBoundaryMakers;
		this.NumberOfHoles = numberOfHoles;
		this.NumberOfRegionalAttrs = numberOfRegionalAttrs;
		this.Segments = segments;
		this.Holes = holes;
	}

	hasSegments() {
		return this.Segments.length > 0;
	}
	count() {
		return this.Vertices.length;
	}
}


class Node {
	constructor(numOfVertices, dimension, numOfAttributes, numOfBoundaryMakers) {
		this.NumOfVertices = numOfVertices;
		this.Dimension = dimension;
		this.NumOfAttributes = numOfAttributes;
		this.NumOfBoundaryMakers = numOfBoundaryMakers;
		this.Vertices = [];
		this.BoundingBox = new BoundingBox();
	}

	addVertex(vertex) {
		this.Vertices.push(vertex);
		this.BoundingBox.update(vertex.X, vertex.Y);
	}
}


class Poly {
	constructor(node, numOfSegments, numOfBoundaryMakers) {
		this.Node = node;
		this.NumOfSegments = numOfSegments;
		this.NumOfBoundaryMakers = numOfBoundaryMakers;
		this.NumberOfHoles = 0;
		this.NumberOfRegionalAttrs = 0;
		this.Segments = [];
		this.Holes = [];
	}

	addSegment(segment) {
		this.Segments.push(segment);
	}

	addHole(hole) {
		this.Holes.push(hole);
	}
}

//A queue used to store encroached subsegments.
//Each subsegment's vertices are stored so that we can check whether a 
//subsegment is still the same.
class BadSubSeg {
	constructor() {
		// An encroached subsegment.
		BadSubSeg.HashSeed = 0;
		this.EncSubSeg;
		// Its two vertices.
		this.SubSegOrg;
		this.SubSegDest;
		//BadSubSeg.HashSeed = 0;
		this.Hash = BadSubSeg.HashSeed++;
	}

	getHash() {
		return this.Hash;
	}
}

// An oriented subsegment.
// Includes a pointer to a subsegment and an orientation. The orientation
// denotes a side of the edge.  There are two possible orientations.
// By convention, the edge is always directed so that the "side" denoted
// is the right side of the edge.
class Osub {
	constructor() {
		this.Segment = null;
		// Ranges from 0 to 1.
		this.Orient = 0;
	}

	clone() {
		var orient = JSON.parse(JSON.stringify(this.Orient));
		var osub = new Osub();
		osub.Orient = orient;
		osub.Segment = this.Segment;
		return osub;
	}

    // Set destination of a subsegment.
	setDest(ptr) {
		this.Segment.Vertices[1 - this.Orient] = ptr;
	}

	// Get the origin of the segment that includes the subsegment.
	segOrg() {
		return this.Segment.Vertices[2 + this.Orient];
	}

    // Get the destination of the segment that includes the subsegment.
	segDest() {
		return this.Segment.Vertices[3 - this.Orient];
	}

    // Set the origin of the segment that includes the subsegment.
	setSegOrg(ptr) {
		this.Segment.Vertices[2 + this.Orient] = ptr;
	}
	// Set the destination of the segment that includes the subsegment.
	setSegDest(ptr) {
		this.Segment.Vertices[3 - this.Orient] = ptr;
	}

	// Reverse the orientation of a subsegment. [sym(ab) -> ba]
	symSelf() {
		this.Orient = 1 - this.Orient;
	}

	// Find adjoining subsegment with the same origin. [pivot(ab) -> a*]
	pivot() {
		return this.Segment.SubSegs[this.Orient];
    }
    
	// Bond two subsegments together. [bond(abc, ba)]
	bond(o2) {
		this.Segment.SubSegs[this.Orient] = o2.clone();
		o2.Segment.SubSegs[o2.Orient] = this.clone();
		return o2;
	}

	// Find next subsegment in sequence. [next(ab) -> b*]
	nextSelf() {
		var self = this.clone();
		self = this.Segment.SubSegs[1 - this.Orient];
	}

	// Set the origin or destination of a subsegment.
	setOrg(ptr) {
		this.Segment.Vertices[this.Orient] = ptr;
	}

	// Set destination of a subsegment.
	setDest(ptr) {
		this.Segment.Vertices[1 - this.Orient] = ptr;
	}

	// Set a subsegment's deallocation.
	static kill(sub) {
		sub.SubSegs[0].Segment = null;
		sub.SubSegs[1].Segment = null;
	}

	triDissolve() {
		this.Segment.Triangles[this.Orient].Triangle = Mesh.dummyTri();
	}
}


// An oriented triangle.
// Includes a pointer to a triangle and orientation.
// The orientation denotes an edge of the triangle. There are
// three possible orientations. By convention, each edge always points
// counterclockwise about the corresponding triangle.
class Otri {
	constructor() {
		this.Plus1Mod3 = [1, 2, 0];
		this.Minus1Mod3 = [2, 0, 1];
		this.Triangle = null;
		// Ranges from 0 to 2.
		this.Orient = 0;
	}
	clone() {
		var plus1Mod3 = JSON.parse(JSON.stringify(this.Plus1Mod3));
		var minus1Mod3 = JSON.parse(JSON.stringify(this.Minus1Mod3));
		var orient = JSON.parse(JSON.stringify(this.Orient));
		var otri = new Otri();
		otri.Plus1Mod3 = plus1Mod3;
		otri.Minus1Mod3 = minus1Mod3;
		otri.Orient = orient;
		otri.Triangle = this.Triangle;
		return otri;
	}

	// Dissolve a bond (from one side).  
	//The other triangle will still think it's connected to 
	// this triangle. The other triangle is being deleted 
	// entirely, or bonded to another triangle, so it doesn't matter.
	dissolve() {
		this.Triangle.Neighbors[this.Orient].Triangle = Mesh.dummyTri();
		this.Triangle.Neighbors[this.Orient].Orient = 0;
	}

	// Copy an oriented triangle.
	copy(o2) {
		o2.Triangle = this.Triangle;
		o2.Orient = this.Orient;
	}

	equal(o2) {
		return (this.Triangle === o2.Triangle && this.Orient === o2.Orient);
	}

	setOrg(vertex) {
		this.Triangle.Vertices[this.Plus1Mod3[this.Orient]] = vertex;
	}

	setDest(vertex) {
		this.Triangle.Vertices[this.Minus1Mod3[this.Orient]] = vertex;
	}
	setApex(vertex) {
		this.Triangle.Vertices[this.Orient] = vertex;
	}

	// Bond two triangles together at the resepective handles. [bond(abc, bad)]
	bond(o2) {
		this.Triangle.Neighbors[this.Orient].Triangle = o2.Triangle;
		this.Triangle.Neighbors[this.Orient].Orient = o2.Orient;
		o2.Triangle.Neighbors[o2.Orient].Triangle = this.Triangle;
		o2.Triangle.Neighbors[o2.Orient].Orient = this.Orient;
		return o2;
	}

	org() {
		return this.Triangle.Vertices[this.Plus1Mod3[this.Orient]];
	}
	dest() {
		return this.Triangle.Vertices[this.Minus1Mod3[this.Orient]];
	}
	apex() {
		return this.Triangle.Vertices[this.Orient];
	}

	// Find the previous edge (clockwise) of a triangle. [lprev(abc) -> cab]
	lPrev(o2) {
		o2.Triangle = this.Triangle;
		o2.Orient = this.Minus1Mod3[this.Orient];
		return o2;
	}

	// Find the next edge clockwise with the same origin. [oprev(abc) -> a*b]
	// oprev() spins clockwise around a vertex; It finds the 
	// next edge with the same origin in the clockwise direction.  This edge is 
	// part of a different triangle.
	oPrev(o2) {
		o2.Triangle = this.Triangle.Neighbors[this.Orient].Triangle;
		o2.Orient = this.Triangle.Neighbors[this.Orient].Orient;

		o2.Orient = this.Plus1Mod3[o2.Orient];
		return o2;
	}

	// Find the previous edge clockwise with the same origin. [oprev(abc) -> a*b]
	oPrevSelf() {
		let tmp = this.Orient;
		this.Orient = this.Triangle.Neighbors[tmp].Orient;
		this.Triangle = this.Triangle.Neighbors[tmp].Triangle;
		this.Orient = this.Plus1Mod3[this.Orient];
	}

	// Find the previous edge (clockwise) of a triangle. [lprev(abc) -> cab]
	lPrevSelf() {
		this.Orient = this.Minus1Mod3[this.Orient];
	}

    // Find the next edge (counterclockwise) of a triangle. [lnext(abc) -> bca]
	lNext(o2) {
		o2.Triangle = this.Triangle;
		o2.Orient = this.Plus1Mod3[this.Orient];
		return o2;
	}

    // Find the next edge counterclockwise with the same origin. [onext(abc) -> ac*]
	// onext() spins counterclockwise around a vertex; It finds 
	// the next edge with the same origin in the counterclockwise direction. This
	// edge is part of a different triangle.
	oNext(o2) {
		o2.Triangle = this.Triangle;
		o2.Orient = this.Minus1Mod3[this.Orient];

		let tmp = o2.Orient;
		o2.Orient = o2.Triangle.Neighbors[tmp].Orient;
		o2.Triangle = o2.Triangle.Neighbors[tmp].Triangle;
		return o2;
	}

	// Find the next edge (counterclockwise) of a triangle. [lnext(abc) -> bca]
	lNextSelf() {
		this.Orient = this.Plus1Mod3[this.Orient];
	}

	// Find the abutting triangle; same edge. [sym(abc) -> ba*]      
	// Note that the edge direction is necessarily reversed, because the handle specified 
	// by an oriented triangle is directed counterclockwise around the triangle.
	sym(o2) {
		o2.Triangle = this.Triangle.Neighbors[this.Orient].Triangle;
		o2.Orient = this.Triangle.Neighbors[this.Orient].Orient;
		return o2;
    }
    
	// Find the abutting triangle; same edge. [sym(abc) -> ba*]
	symSelf() {
		let tmpOrient = this.Orient;
		this.Orient = this.Triangle.Neighbors[tmpOrient].Orient;
		this.Triangle = this.Triangle.Neighbors[tmpOrient].Triangle;
	}

	//Find the next edge counterclockwise with the same origin. [onext(abc) -> ac*]
	oNextSelf() {
		this.Orient = this.Minus1Mod3[this.Orient];
		this.symSelf();
	}

	//Find the next edge clockwise with the same destination. [dprev(abc) -> cb*]
	dPrevSelf() {
		this.lNextSelf();
		this.symSelf();
	}

	// Set a triangle's deallocation.
	static kill(tria) {
		tria.Neighbors[0].Triangle = null;
		tria.Neighbors[2].Triangle = null;
    }
    
	// Finds a subsegment abutting a triangle.
	segPivot() {
		return this.Triangle.SubSegs[this.Orient];
    }
    
	// Bond a triangle to a subsegment.    
	segBond(os) {
		this.Triangle.SubSegs[this.Orient] = os.clone();
		os.Segment.Triangles[os.Orient] = this;
		return os;
	}

	// Dissolve a bond (from the triangle side).
	segDissolve() {
		this.Triangle.SubSegs[this.Orient].Segment = Mesh.dummySub();
	}

	//Check a triangle's deallocation.
	static isDead(tria) {
		return tria.Neighbors[0].Triangle == null;
	}

	// Test a triangle for viral infection.
	isInfected() {
		return this.Triangle.Infected;
	}

	// Infect a triangle with the virus.
	infect() {
		this.Triangle.Infected = true;
	}

	// Cure a triangle from the virus.
	uninfect() {
		this.Triangle.Infected = false;
	}
}

// The subsegment data structure.
// Each subsegment contains two pointers to adjoining subsegments, plus
// four pointers to vertices, plus two pointers to adjoining triangles,
// plus one boundary marker.
class Segment {

	constructor(segmentNum, endPoint1, endPoint2, boundaryMarker) {
		this.SegmentNum = segmentNum;
		this.EndPoint1 = endPoint1;
		this.EndPoint2 = endPoint2;
		this.BoundaryMarker = boundaryMarker;
		this.SubSegs = [2];
		this.SubSegs[0] = new Osub();
		this.SubSegs[1] = new Osub();
		this.SubSegs[0].Segment = Mesh.dummySub();
		this.SubSegs[1].Segment = Mesh.dummySub();
		this.Vertices = [4];
		this.Vertices[0] = new Vertex();
		this.Vertices[1] = new Vertex();
		this.Vertices[2] = new Vertex();
		this.Vertices[3] = new Vertex();
		this.Triangles = [2];
		this.Triangles[0] = new Otri();
		this.Triangles[1] = new Otri();
		this.Triangles[0].Triangle = Mesh.dummyTri();
		this.Triangles[1].Triangle = Mesh.dummyTri();
		this.BoundaryMarker = 0;
		this.Hash = 0;
	}

	p0() {
		return this.Vertices[0].Id;
	}

	p1() {
		return this.Vertices[1].Id;
	}
}

// The triangle data structure.
// Each triangle contains three pointers to adjoining triangles, plus three 
// pointers to vertices, plus three pointers to subsegments (declared below;
// these pointers are usually 'dummysub'). It may or may not also contain 
// user-defined attributes and/or a floating-point "area constraint".
class Triangle {
	constructor() {
		this.Id;
		this.Neighbors = [3];
		this.Vertices = [3];
		this.SubSegs = [3];
		this.Region = 0;
		this.Area = 0;
		this.Infected;
		this.Vertices[0] = new Vertex();
		this.Vertices[1] = new Vertex();
		this.Vertices[2] = new Vertex();
		this.Neighbors[0] = new Otri();
		this.Neighbors[1] = new Otri();
		this.Neighbors[2] = new Otri();
		this.Neighbors[0].Triangle = Mesh.dummyTri();
		this.Neighbors[1].Triangle = Mesh.dummyTri();
		this.Neighbors[2].Triangle = Mesh.dummyTri();
		this.SubSegs[0] = new Osub();
		this.SubSegs[1] = new Osub();
		this.SubSegs[2] = new Osub();
		this.SubSegs[0].Segment = Mesh.dummySub();
		this.SubSegs[1].Segment = Mesh.dummySub();
		this.SubSegs[2].Segment = Mesh.dummySub();
	}

	p0() {
		if (this.Vertices[0] === null) {
			return -1;
		} else {
			return this.Vertices[0].Id;
		}
	}

	p1() {
		if (this.Vertices[1] === null) {
			return -1;
		} else {
			return this.Vertices[1].Id;
		}
	}
	p2() {
		if (this.Vertices[2] === null) {
			return -1;
		} else {
			return this.Vertices[2].Id;
		}
	}

	n0() {
		return this.Neighbors[0].Triangle.Id;
	}

	n1() {
		return this.Neighbors[1].Triangle.Id;
	}

	n2() {
		return this.Neighbors[2].Triangle.Id;
	}
}

class Primitives {

	constructor() {
		Primitives.Splitter = 0;
		Primitives.Epsilon = 0;
		Primitives.CCWErrBoundA = 0;
		Primitives.ICCErrBoundA = 0;
	}

	static exactInit() {
		let half = 0.5;
		let check = 1.0;
		let lastCheck = 0;
		let every_other = true;
		this.Splitter = 1.0;
		this.Epsilon = 1.0;

        do {
			lastCheck = check;
			this.Epsilon *= half;
			if (every_other === true) {
				this.Splitter *= 2.0;
			}
			every_other = !every_other;
			check = 1.0 + this.Epsilon;
		} while ((check !== 1.0) && (check !== lastCheck));
		this.Splitter += 1.0;
		// Error bounds for orientation and incircle tests.
		this.CCWErrBoundA = (3.0 + 16.0 * this.Epsilon) * this.Epsilon;
		this.ICCErrBoundA = (10.0 + 96.0 * this.Epsilon) * this.Epsilon;
	}

	//Returns twice the area of the oriented triangle (a,b,c), i.e., the
	//area is positive if the triangle is oriented counterclockwise
	static triArea(a, b, c) {
		return (b.X - a.X) * (c.Y - a.Y) - (b.Y - a.Y) * (c.X - a.X);
	}

    // Check if the point pd lies inside the circle passing through pa, pb, and pc. The 
	// points pa, pb, and pc must be in counterclockwise order, or the sign of the result 
	// will be reversed.
	// Return a positive value if the point pd lies inside the circle passing through 
	// pa, pb, and pc; a negative value if it lies outside; and zero if the four points 
	// are cocircular.
	// Uses exact arithmetic if necessary to ensure a correct answer.  The
	// result returned is the determinant of a matrix.  This determinant is
	// computed adaptively, in the sense that exact arithmetic is used only to
	// the degree it is needed to ensure that the returned value has the
	// correct sign.
	static inCircle(pa, pb, pc, pd) {
		let adx = 0;
		let bdx = 0;
		let cdx = 0;
		let ady = 0;
		let bdy = 0;
		let cdy = 0;
		let bdxcdy = 0;
		let cdxbdy = 0;
		let cdxady = 0;
		let adxcdy = 0;
		let adxbdy = 0;
		let bdxady = 0;
		let alift = 0;
		let blift = 0;
		let clift = 0;
		let det = 0;
		let permanent = 0;
		let errbound = 0;

		adx = pa.X - pd.X;
		bdx = pb.X - pd.X;
		cdx = pc.X - pd.X;
		ady = pa.Y - pd.Y;
		bdy = pb.Y - pd.Y;
		cdy = pc.Y - pd.Y;

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;
		alift = adx * adx + ady * ady;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;
		blift = bdx * bdx + bdy * bdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;
		clift = cdx * cdx + cdy * cdy;

		det = alift * (bdxcdy - cdxbdy) +
			blift * (cdxady - adxcdy) +
			clift * (adxbdy - bdxady);

		permanent = (Math.abs(bdxcdy) + Math.abs(cdxbdy)) * alift +
			(Math.abs(cdxady) + Math.abs(adxcdy)) * blift +
			(Math.abs(adxbdy) + Math.abs(bdxady)) * clift;
		errbound = this.ICCErrBoundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}
		return Primitives.inCircleDecimal(pa, pb, pc, pd);
	}

	static inCircleDecimal(pa, pb, pc, pd) {
		let adx = 0;
		let bdx = 0;
		let cdx = 0;
		let ady = 0;
		let bdy = 0;
		let cdy = 0;
		let bdxcdy = 0;
		let cdxbdy = 0;
		let cdxady = 0;
		let adxcdy = 0;
		let adxbdy = 0;
		let bdxady = 0;
		let alift = 0;
		let blift = 0;
		let clift = 0;

		adx = pa.X - pd.X;
		bdx = pb.X - pd.X;
		cdx = pc.X - pd.X;
		ady = pa.Y - pd.Y;
		bdy = pb.Y - pd.Y;
		cdy = pc.Y - pd.Y;

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;
		alift = adx * adx + ady * ady;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;
		blift = bdx * bdx + bdy * bdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;
		clift = cdx * cdx + cdy * cdy;

		return alift * (bdxcdy - cdxbdy) +
			blift * (cdxady - adxcdy) +
			clift * (adxbdy - bdxady);

	}


	// Check, if the three points appear in counterclockwise order. The result is 
	// also a rough approximation of twice the signed area of the triangle defined 
	// by the three points.
	// Return a positive value if the points pa, pb, and pc occur in 
	// counterclockwise order; a negative value if they occur in clockwise order; 
	// and zero if they are collinear
	static ccw(pa, pb, pc) {
		let detLeft = (pa.X - pc.X) * (pb.Y - pc.Y);
		let detRight = (pa.Y - pc.Y) * (pb.X - pc.X);
		let det = detLeft - detRight;
		return det;
	}

	static rightOf(x, triangle) {
		if (triangle === null) {
			return false;
		}
		return this.ccw(x, triangle.dest(), triangle.org());
	}

}
class TriangleLocator {
	constructor(mesh) {
		this.Mesh = mesh;
		// Pointer to a recently visited triangle. Improves point location if
		// proximate vertices are inserted sequentially.
		this.RecentTri = new Otri();
		this.Sampler = new Sampler();
	}

	update(otri) {
		otri.copy(this.RecentTri);
	}

	reset() {
		this.RecentTri.Triangle = null; // No triangle has been visited yet.
	}

	locate(searchPoint, searchTri) {
		let sampleTri = new Otri();
		let tOrg = null;
		let tDest = null;

		let searchDist = 0;
		let dist = 0;
		let ahead = 0;
		// Record the distance from the suggested starting triangle to the
		// point we seek.
		tOrg = searchTri.org();
		searchDist = (searchPoint.X - tOrg.X) * (searchPoint.X - tOrg.X) + (searchPoint.Y - tOrg.Y) * (searchPoint.Y - tOrg.Y);
		// If a recently encountered triangle has been recorded and has not been
		// deallocated, test it as a good starting point.        
		if (this.RecentTri.Triangle != null) {
			if (!Otri.isDead(this.RecentTri.Triangle)) {
				tOrg = this.RecentTri.org();
				if (tOrg.X == searchPoint.X && tOrg.Y == searchPoint.Y) {
					this.RecentTri.copy(searchTri);
					return LocateResult.OnVertex;
				}
				dist = (searchPoint.X - tOrg.X) * (searchPoint.X - tOrg.X) + (searchPoint.Y - tOrg.Y) * (searchPoint.Y - tOrg.Y);
				if (dist < searchDist) {
					this.RecentTri.copy(searchTri);
					searchDist = dist;
				}
			}
		}

		this.Sampler.update(this.Mesh);
		var samples = this.Sampler.getSamples(this.Mesh);
		for (var key of samples) {
			sampleTri.Triangle = this.Mesh.Triangles.get(key);
			if (!Otri.isDead(sampleTri.Triangle)) {
				tOrg = sampleTri.org();
				dist = (searchPoint.X - tOrg.X) * (searchPoint.X - tOrg.X) + (searchPoint.Y - tOrg.Y) * (searchPoint.Y - tOrg.Y);
			}
			if (dist < searchDist) {
				sampleTri.copy(searchTri);
				searchDist = dist;
			}
		}

		tOrg = searchTri.org();
		tDest = searchTri.dest();
		// Check the starting triangle's vertices.
		if (tOrg.X === searchPoint.X && tOrg.Y === searchPoint.Y) {
			return LocateResult.OnVertex;
		}
		if (tDest.X === searchPoint.X && tDest.Y === searchPoint.Y) {
			searchTri.lNextSelf();
			return LocateResult.OnVertex;
		}
		// Orient 'searchtri' to fit the preconditions of calling preciselocate()
		ahead = Primitives.ccw(tOrg, tDest, searchPoint);
		if (ahead < 0.0) {
			// Turn around so that 'searchpoint' is to the left of the
			// edge specified by 'searchtri'.
			searchTri.symSelf();
		} else if (ahead === 0.0) {
			// Check if 'searchpoint' is between 'torg' and 'tdest'.
			if (((tOrg.X < searchPoint.X) === (searchPoint.X < tDest.X)) &&
				((tOrg.Y < searchPoint.Y) === (searchPoint.Y < tDest.Y))) {
				return LocateResult.OnEdge;
			}
		}
		return this.preciseLocate(searchPoint, searchTri, false);
	}

	preciseLocate(searchPoint, searchTri, stopAtSubsegment) {
		let backTrackTri = new Otri();
		let checkEdge = new Osub();
		let destOrient = 0;
		let orgOrient = 0;
		let fOrg = searchTri.org();
		let fDest = searchTri.dest();
		let fApex = searchTri.apex();
		let moveLeft = false;
		while (true) {
			// Check whether the apex is the point we seek.
			if (fApex.X === searchPoint.X && fApex.Y === searchPoint.Y) {
				searchTri.lPrevSelf();
				return LocateResult.OnVertex;
			}
			// Does the point lie on the other side of the line defined by the
			// triangle edge opposite the triangle's destination?
			destOrient = Primitives.ccw(fOrg, fApex, searchPoint);
			// Does the point lie on the other side of the line defined by the
			// triangle edge opposite the triangle's origin?
			orgOrient = Primitives.ccw(fApex, fDest, searchPoint);
			if (destOrient > 0.0) {
				if (orgOrient > 0.0) {
					// Move left if the inner product of (fapex - searchpoint) and
					// (fdest - forg) is positive.  This is equivalent to drawing
					// a line perpendicular to the line (forg, fdest) and passing
					// through 'fapex', and determining which side of this line
					// 'searchpoint' falls on.
					moveLeft = (fApex.X - searchPoint.X) * (fDest.X - fOrg.X) +
						(fApex.Y - searchPoint.Y) * (fDest.Y - fOrg.Y) > 0.0;
				} else {
					moveLeft = true;
				}
			} else {
				if (orgOrient > 0.0) {
					moveLeft = false;
				} else {
					// The point we seek must be on the boundary of or inside this
					// triangle.
					if (destOrient === 0.0) {
						searchTri.lPrevSelf();
						return LocateResult.OnEdge;
					}
					if (orgOrient === 0.0) {
						searchTri.lNextSelf();
						return LocateResult.OnEdge;
					}
					return LocateResult.InTriangle;
				}
			}
			// Move to another triangle. Leave a trace 'backTrackTri' in case
			// floating-point roundoff or some such bogey causes us to walk
			// off a boundary of the triangulation.
			if (moveLeft === true) {
				backTrackTri = searchTri.lPrev(backTrackTri);
				fDest = fApex;
			} else {
				backTrackTri = searchTri.lNext(backTrackTri);
				fOrg = fApex;
			}
			searchTri = backTrackTri.sym(searchTri);
			if (this.Mesh.CheckSegments && stopAtSubsegment) {
				// Check for walking through a subsegment.
				checkEdge = backTrackTri.segPivot();
				if (checkEdge.Segment !== Mesh.dummySub()) {
					// Go back to the last triangle.
					backTrackTri.copy(searchTri);
					return LocateResult.Outside;
				}
			}
			//Check for walking right out of the triangulation.
			if (searchTri.Triangle === Mesh.dummyTri()) {
				// Go back to the last triangle.
				backTrackTri.copy(searchTri);
				return LocateResult.Outside;
			}
			fApex = searchTri.apex();
		}
	}
}

class Mesh {
	constructor(inputGeometry) {
		// The 'triangle' that occupies all of 'outer space'.
		Mesh._DummyTri = null;
		// The omnipresent subsegment. Referenced by any triangle or subsegment
		// that isn't really connected to a subsegment at that location.
		Mesh._DummySub = null;
		this.Vertices = new Map();
		this.Triangles = new Map();
		this.SubSegs = new Map();
		this.Holes = [];
		this.Regions = [];
		this.Bounds = inputGeometry.BoundingBox;

		this.InVertices; // Number of input vertices.
		this.InElements; // Number of input triangles.
		this.InSegments; // Number of input segments.
		this.Undeads; // Number of input vertices that don't appear in the mesh.
		this.Edges; // Number of output edges.
		this.Mesh_Dim; // Dimension 
		this.NExtras; // Number of attributes per vertex.
		this.HullSize; // Number of edges in convex hull.
		this.SteinerLeft;
		this.CheckSegments; // Are there segments in the triangulation yet?


		// Triangular bounding box vertices.
		this.InfVertex1 = null;
		this.InfVertex2 = null;
		this.InfVertex3 = null;
		this.CheckSegments = false;
		this.Locator = new TriangleLocator(this);
		// Number of edges in convex hull.
		this.HullSize = 0;
		// Hash seeds (should belong to mesh instance)
		this.Hash_Vtx = 0;
		this.Hash_Seg = 0;
		this.Hash_Tri = 0;
		this.Behavior = new Behavior();
		this.Numbering;

		Primitives.exactInit();
		if (Mesh.dummyTri() === null) {
			// Initialize static dummy triangle and subseg.
			this.dummyInit();
			for (var i = 0; i < inputGeometry.Vertices.length; i++) {
				this.Vertices.set(inputGeometry.Vertices[i].Id, inputGeometry.Vertices[i]);
			}
		}
	}

	static dummyTri() {
		return Mesh._DummyTri;
	}
	static setDummyTri(value) {
		Mesh._DummyTri = value;
	}

	static dummySub() {
		return Mesh._DummySub;
	}
	static setDummySub(value) {
		Mesh._DummySub = value;
	}

	delaunay() {
		let hullEdges = 0;
		let incremental = new Incremental();
		incremental.triangulate(this);
		// The input vertices may all be collinear, so there are 
		// no triangles.
		return (this.Triangles.length === 0) ? 0 : hullEdges;
	}
	resetData() {
		this.Vertices.clear();
		this.Triangles.clear();
		this.SubSegs.clear();
		this.Holes = [];
		this.Regions = [];

		this.Hash_Vtx = 0;
		this.Hash_Seg = 0;
		this.Hash_Tri = 0;

		this.HullSize = 0;
		this.Edges = 0;

		this.reset();
		this.Locator.reset();
	}

	reset() {
		this.Undeads = 0; // No eliminated input vertices yet.
		this.CheckSegments = false; // There are no segments in the triangulation yet.
    }
    
	triangulate(input) {
		this.resetData();
		this.Behavior.Poly = input.hasSegments();
		if (this.Behavior.Poly === false) {
			this.Behavior.VarArea = false;
		}

		this.Behavior.UseRegions = false;
		this.SteinerLeft = this.Behavior.SteinerPoints;
		this.transferNodes(input);
		this.HullSize = this.delaunay(); // Triangulate the vertices.

		// Ensure that no vertex can be mistaken for a triangular bounding
		// box vertex in insertvertex().
		this.InfVertex1 = null;
		this.InfVertex2 = null;
		this.InfVertex3 = null;

		if (this.Behavior.UseSegments === true) {
			this.CheckSegments = true;

			// Insert PSLG segments and/or convex hull segments.
			this.formSkeleton(input);
		}
		if (this.Behavior.Poly === true && this.Triangles.size > 0) {
			// Copy holes
			for (var item of input.Holes) {
				this.Holes.push(item);
			}

            // Carve out holes and concavities.
			let c = new Carver(this);
			c.carveHoles();
		} else {
			// Without a PSLG, there can be no holes or regional attributes
			// or area constraints. The following are set to zero to avoid
			// an accidental free() later.
			//
			this.Holes.clear();
			this.Regions.clear();
		}
		// Calculate the number of edges.
		this.Edges = (3 * this.Triangles.size + this.HullSize) / 2;
	}

	formSkeleton(input) {
		let endPoint1 = null;
		let endPoint2 = null;
		let end1 = 0;
		let end2 = 0;
		let boundMarker = 0;
		this.InSegments = 0;

		if (this.Behavior.Poly === true) {
			// If the input vertices are collinear, there is no triangulation,
			// so don't try to insert segments.
			if (this.Triangles.length === 0) {
				return;
			}

			// If segments are to be inserted, compute a mapping
			// from vertices to triangles.
			if (input.hasSegments() === true) {
				this.makeVertexMap();
			}
			boundMarker = 0;
			// Read and insert the segments.
			for (var seg of input.Segments) {
				this.InSegments++;
				end1 = seg.P0;
				end2 = seg.P1;
				boundMarker = seg.Boundary;
				if (end1 < 0 || end1 >= this.InVertices) {

				} else if (end2 < 0 || end2 >= this.InVertices) {

				} else {
					// Find the vertices numbered 'end1' and 'end2'.
					endPoint1 = this.Vertices.get(end1);
					endPoint2 = this.Vertices.get(end2);
					if (endPoint1.X === endPoint2.X && endPoint1.Y === endPoint2.Y) {

					} else {
						this.insertSegment(endPoint1, endPoint2, boundMarker);
					}
				}
			}
		}
		if (this.Behavior.Convex === true || this.Behavior.Poly === false) {
			// Enclose the convex hull with subsegments.
			this.markHull();
		}
	}

	// Insert a PSLG segment into a triangulation.
	insertSegment(endPoint1, endPoint2, newMark) {
		let searchTri1 = new Otri();
		let searchTri2 = new Otri();
		let checkVertex = null;
		// Find a triangle whose origin is the segment's first endpoint.
		searchTri1 = endPoint1.Tri;
		if (searchTri1.Triangle !== null) {
			checkVertex = searchTri1.org();
		}
		if (!checkVertex.equals(endPoint1)) {
			// Find a boundary triangle to search from.
			searchTri1.Triangle = Mesh.dummyTri();
			searchTri1.Orient = 0;
			searchTri1.symSelf();
			// Search for the segment's first endpoint by point location.
			if (this.Locator.locate(endPoint1, searchTri1) !== LocateResult.OnVertex) {
				alert("Unable to locate PSLG vertex in triangulation");
			}
		}
		// Remember this triangle to improve subsequent point location.
		this.Locator.update(searchTri1);
		// Scout the beginnings of a path from the first endpoint
		// toward the second.
		if (this.scoutSegment(searchTri1, endPoint2, newMark) === true) {
			// The segment was easily inserted.
			return;
		}
		// The first endpoint may have changed if a collision with an intervening
		// vertex on the segment occurred.
		endPoint1 = searchTri1.org();
		// Find a triangle whose origin is the segment's second endpoint.
		checkVertex = null;
		searchTri2 = endPoint2.Tri;
		if (searchTri2.Triangle !== null) {
			checkVertex = searchTri2.org();
		}
		if (checkVertex !== endPoint2) {
			// Find a boundary triangle to search from.
			searchTri2.Triangle = Mesh.dummyTri();
			searchTri2.Orient = 0;
			searchTri2.symSelf();
			if (this.Locator.locate(endPoint2, searchTri2) !== LocateResult.OnVertex) {
				alert('Unable to locate PSLG vertex in triangulation.');
			}
		}
		// Remember this triangle to improve subsequent point location.
		this.Locator.update(searchTri2);
		// Scout the beginnings of a path from the second endpoint
		// toward the first.
		if (this.scoutSegment(searchTri2, endPoint1, newMark)) {
			// The segment was easily inserted.
			return;
		}
		// The second endpoint may have changed if a collision with an intervening
		// vertex on the segment occurred.
		endPoint2 = searchTri2.org();
		// Insert the segment directly into the triangulation.
		this.constrainedEdge(searchTri1, endPoint2, newMark);
	}

	// Cover the convex hull of a triangulation with subsegments.
	markHull() {
		let hullTri = new Otri();
		let nextTri = new Otri();
		let startTri = new Otri();

		// Find a triangle handle on the hull.
		hullTri.Triangle = Mesh.dummyTri();
		hullTri.Orient = 0;
		hullTri.symSelf();
		// Remember where we started so we know when to stop.
		hullTri.copy(startTri);
		// Go once counterclockwise around the convex hull.
		do {
			// Create a subsegment if there isn't already one here.
			this.insertSubSeg(hullTri, 1);
			// To find the next hull edge, go clockwise around the next vertex.
			hullTri.lNextSelf();
			nextTri = hullTri.oPrev(nextTri);
			while (nextTri.Triangle !== Mesh.dummyTri()) {
				nextTri.copy(hullTri);
				nextTri = hullTri.oPrev(nextTri);
			}
		} while (hullTri.equal(startTri) === false)
	}

	// Force a segment into a constrained Delaunay triangulation by deleting the 
	// triangles it intersects, and triangulating the polygons that form on each 
	// side of it.
	// Generates a single subsegment connecting 'endpoint1' to 'endpoint2'.
	// The triangle 'starttri' has 'endpoint1' as its origin.  'newmark' is the
	// boundary marker of the segment.
	//
	// To insert a segment, every triangle whose interior intersects the
	// segment is deleted. The union of these deleted triangles is a polygon
	// (which is not necessarily monotone, but is close enough), which is
	// divided into two polygons by the new segment. This routine's task is
	// to generate the Delaunay triangulation of these two polygons.
	//
	// This routine's behavior ia a two-step process.  The
	// first step is to walk from endpoint1 to endpoint2, flipping each edge
	// encountered.  This step creates a fan of edges connected to endpoint1,
	// including the desired edge to endpoint2. The second step enforces the
	// Delaunay condition on each side of the segment in an incremental manner:
	// proceeding along the polygon from endpoint1 to endpoint2 (this is done
	// independently on each side of the segment), each vertex is "enforced"
	// as if it had just been inserted, but affecting only the previous
	// vertices. The result is the same as if the vertices had been inserted
	// in the order they appear on the polygon, so the result is Delaunay.
	//
	// ConstrainedEdge() interleaves these two steps. The procedure
	// walks from endpoint1 to endpoint2, and each time an edge is encountered
	// and flipped, the newly exposed vertex (at the far end of the flipped
	// edge) is "enforced" upon the previously flipped edges, usually affecting
	// only one side of the polygon (depending upon which side of the segment
	// the vertex falls on).
	//
	// Although the polygon is not necessarily monotone, it can be
	// triangulated in a manner similar to the stack-based algorithms for
	// monotone polygons. For each reflex vertex (local concavity) of the
	// polygon, there will be an inverted triangle formed by one of the edge
	// flips. (An inverted triangle is one with negative area - that is, its
	// vertices are arranged in clockwise order - and is best thought of as a
	// wrinkle in the fabric of the mesh.)  Each inverted triangle can be
	// thought of as a reflex vertex pushed on the stack, waiting to be fixed
	// later.
	//
	// A reflex vertex is popped from the stack when a vertex is inserted that
	// is visible to the reflex vertex. (However, if the vertex behind the
	// reflex vertex is not visible to the reflex vertex, a new inverted
	// triangle will take its place on the stack.) These details are handled
	// by the DelaunayFixup() routine above.
	constrainedEdge(startTri, endPoint2, newMark) {
		let fixUpTri = new Otri();
		let fixUpTri2 = new Otri();
		let crossSubSeg = new Osub();
		let endPoint1 = null;
		let farVertex = null;
		let area = 0;
		let collision = false;
		let done = false;

		endPoint1 = startTri.org();
		fixUpTri = startTri.lNext(fixUpTri);
		this.flip(fixUpTri);
		// 'collision' indicates whether we have found a vertex directly
		// between endpoint1 and endpoint2.
		do {
			farVertex = fixUpTri.org();
			// 'farvertex' is the extreme point of the polygon we are "digging"
			//  to get from endpoint1 to endpoint2.
			if (farVertex.X === endPoint2.X && farVertex.Y === endPoint2.Y) {
				fixUpTri2 = fixUpTri.oPrev(fixUpTri2);
				// Enforce the Delaunay condition around endpoint2.
				this.delaunayFixup(fixUpTri, false);
				this.delaunayFixup(fixUpTri2, true);
				done = true;
			} else {
				// Check whether farvertex is to the left or right of the segment being
				// inserted, to decide which edge of fixuptri to dig through next.
				area = Primitives.ccw(endPoint1, endPoint2, farVertex);
				if (area === 0.0) {
					// We've collided with a vertex between endpoint1 and endpoint2.
					collision = true;
					fixUpTri2 = fixUpTri.oPrev(fixUpTri2);
					// Enforce the Delaunay condition around farvertex.
					this.delaunayFixup(fixUpTri, false);
					this.delaunayFixup(fixUpTri2, false);
					done = true;
				} else {
					if (area > 0.0) {
						// farvertex is to the left of the segment.
						fixUpTri2 = fixUpTri.oPrev(fixUpTri2);
						// Enforce the Delaunay condition around farvertex, on the
						// left side of the segment only.
						this.delaunayFixup(fixUpTri2, true);
						// Flip the edge that crosses the segment. After the edge is
						// flipped, one of its endpoints is the fan vertex, and the
						// destination of fixuptri is the fan vertex.
						fixUpTri.lPrevSelf();
					} else {
						// farvertex is to the right of the segment.
						this.delaunayFixup(fixUpTri, false);
						// Flip the edge that crosses the segment. After the edge is
						// flipped, one of its endpoints is the fan vertex, and the
						// destination of fixuptri is the fan vertex.
						fixUpTri.oPrevSelf();
					}
					// Check for two intersecting segments.
					crossSubSeg = fixUpTri.segPivot();
					if (crossSubSeg.Segment === Mesh.dummySub()) {
						this.flip(fixUpTri); // May create inverted triangle at left.
					} else {
						// We've collided with a segment between endpoint1 and endpoint2.
						collision = true;
						// Insert a vertex at the intersection.
						this.segmentIntersection(fixUpTri, crossSubSeg, endPoint2);
						done = true;
					}
				}
			}
		} while (done === false);
		// Insert a subsegment to make the segment permanent.
		this.insertSubSeg(fixUpTri, newMark);
		// If there was a collision with an interceding vertex, install another
		// segment connecting that vertex with endpoint2.
		if (collision === true) {
			// Insert the remainder of the segment.
			if (scoutSegment(fixUpTri, endPoint2, newMark) === false) {
				this.constrainedEdge(fixUpTri, endPoint2, newMark);
			}
		}
	}

	// Transform two triangles to two different triangles by flipping an edge 
	// counterclockwise within a quadrilateral.
	// The original triangles, abc and bad, are oriented so that the
	// shared edge ab lies in a horizontal plane, with the vertex b on the left
	// and the vertex a on the right. The vertex c lies below the edge, and
	// the vertex d lies above the edge. The 'flipedge' handle holds the edge
	// ab of triangle abc, and is directed left, from vertex a to vertex b.
	//
	// The triangles abc and bad are deleted and replaced by the triangles cdb
	// and dca.  The triangles that represent abc and bad are NOT deallocated;
	// they are reused for dca and cdb, respectively.  Hence, any handles that
	// may have held the original triangles are still valid, although not
	// directed as they were before.
	//
	// Upon completion of this routine, the 'flipedge' handle holds the edge
	// dc of triangle dca, and is directed down, from vertex d to vertex c.
	// (The two triangles have rotated counterclockwise.)
	//
	// This transformation is geometrically valid only if the
	// quadrilateral adbc is convex.  This transformation is
	// valid only if there is not a subsegment between the triangles abc and
	// bad.  This routine does not check either of these preconditions, and
	// it is the responsibility of the calling routine to ensure that they are
	// met.
	// 
	// A "local transformation" replaces a small set of triangles with another
	// set of triangles.  This may or may not involve inserting or deleting a
	// vertex.
	//
	// The term "casing" is used to describe the set of triangles that are
	// attached to the triangles being transformed, but are not transformed
	// themselves.

	flip(flipEdge) {
		let botLeft = new Otri();
		let botRight = new Otri();
		let topLeft = new Otri();
		let topRight = new Otri();
		let top = new Otri();
		let botLCasing = new Otri();
		let botRCasing = new Otri();
		let topLCasing = new Otri();
		let topRCasing = new Otri();
		let botLSubSeg = new Otri();
		let botRSubSeg = new Otri();
		let topLSubSeg = new Otri();
		let topRSubSeg = new Otri();
		let leftVertex = null;
		let rightVertex = null;
		let botVertex = null;
		let farVertex = null;

		// Identify the vertices of the quadrilateral.
		rightVertex = flipEdge.org();
		leftVertex = flipEdge.dest();
		botVertex = flipEdge.apex();
		top = flipEdge.sym(top);

		farVertex = top.apex();

		// Identify the casing of the quadrilateral.
		topLeft = top.lPrev(topLeft);
		topLCasing = topLeft.sym(topLCasing);
		topRight = top.lNext(topRight);
		topRCasing = topRight.sym(topRCasing);
		botLeft = flipEdge.lNext(botLeft);
		botLCasing = botLeft.sym(botLCasing);
		botRight = flipEdge.lPrev(botRight);
		botRCasing = botRight.sym(botRCasing);

		// Rotate the quadrilateral one-quarter turn counterclockwise.
		botLCasing = topLeft.bond(botLCasing);
		botRCasing = botLeft.bond(botRCasing);
		topRCasing = botRight.bond(topRCasing);
		topLCasing = topRight.bond(topLCasing);

		if (this.CheckSegments === true) {
			// Check for subsegments and rebond them to the quadrilateral.
			topLSubSeg = topLeft.segPivot();
			botLSubSeg = botLeft.segPivot();
			botRSubSeg = botRight.segPivot();
			topRSubSeg = topRight.segPivot();

			if (topLSubSeg.Segment === Mesh.dummySub()) {
				topRight.segDissolve();
			} else {
				topLSubSeg = topRight.segBond(topLSubSeg);
			}

			if (botLSubSeg.Segment === Mesh.dummySub()) {
				topLeft.segDissolve();
			} else {
				botLSubSeg = topLeft.segBond(botLSubSeg);
			}

			if (botRSubSeg.Segment === Mesh.dummySub()) {
				botLeft.segDissolve();
			} else {
				botRSubSeg = botLeft.segBond(botRSubSeg);
			}

			if (topRSubSeg.Segment === Mesh.dummySub()) {
				botRight.segDissolve();
			} else {
				topRSubSeg = botRight.segBond(topRSubSeg);
			}
		}

		// New vertex assignments for the rotated quadrilateral.
		flipEdge.setOrg(farVertex);
		flipEdge.setDest(botVertex);
		flipEdge.setApex(rightVertex);
		top.setOrg(botVertex);
		top.setDest(farVertex);
		top.setApex(leftVertex);
	}

	// Scout the first triangle on the path from one endpoint to another, and check 
	// for completion (reaching the second endpoint), a collinear vertex, or the 
	// intersection of two segments.
	// Returns true if the entire segment is successfully inserted, and false 
	// if the job must be finished by ConstrainedEdge().
	// If the first triangle on the path has the second endpoint as its
	// destination or apex, a subsegment is inserted.
	//
	// If the first triangle on the path has a destination or apex that lies on
	// the segment, a subsegment is inserted connecting the first endpoint to
	// the collinear vertex, and the search is continued from the collinear
	// vertex.
	//
	// If the first triangle on the path has a subsegment opposite its origin,
	// then there is a segment that intersects the segment being inserted.
	// Their intersection vertex is inserted, splitting the subsegment.
	scoutSegment(searchTri, endPoint2, newMark) {
		let crossTri = new Otri();
		let crossSubSeg = new Osub();
		let leftVertex = null;
		let rightVertex = null;
		let collinear = null;

		collinear = this.findDirection(searchTri, endPoint2);
		rightVertex = searchTri.dest();
		leftVertex = searchTri.apex();

		if ((leftVertex.X === endPoint2.X && leftVertex.Y === endPoint2.Y) ||
			(rightVertex.X === endPoint2.X && rightVertex.Y === endPoint2.Y)) {
			// The segment is already an edge in the mesh.
			if (leftVertex.X === endPoint2.X && leftVertex.Y === endPoint2.Y) {
				searchTri.lPrevSelf();
			}
			this.insertSubSeg(searchTri, newMark);
			return true;
		} else if (collinear === FindDirectionResult.Leftcollinear) {
			// We've collided with a vertex between the segment's endpoints.
			// Make the collinear vertex be the triangle's origin.
			searchTri.lPrevSelf();
			this.insertSegment(searchTri, endPoint2, newMark);
			return this.scoutSegment(searchTri, endPoint2, newMark);
		} else if (collinear === FindDirectionResult.Rightcollinear) {
			// We've collided with a vertex between the segment's endpoints.
			this.insertSegment(searchTri, newMark);
			// Make the collinear vertex be the triangle's origin.
			searchTri.lNextSelf();
			// Insert the remainder of the segment.
			return this.scoutSegment(searchTri, endPoint2, newMark);
		} else {
			searchTri.lNext(crossTri);
			crossSubSeg = crossTri.segPivot();
			// Check for a crossing segment.
			if (crossSubSeg.Segment === Mesh.dummySub()) {
				return false;
			} else {
				// Insert a vertex at the intersection.
				this.segmentIntersection(crossTri, crossSubSeg, endPoint2);
				crossTri.copy(searchTri);
				this.insertSubSeg(searchTri, newMark);
				// Insert the remainder of the segment
				return this.scoutSegment(searchTri, endPoint2, newMark);
			}

		}
	}

	// Enforce the Delaunay condition at an edge, fanning out recursively from 
	// an existing vertex.
	// leftside indicates whether or not fixuptri is to the left of 
	// the segment being inserted. 
	// This is a support routine for inserting segments into a constrained
	// Delaunay triangulation.
	//
	// The origin of fixuptri is treated as if it has just been inserted, and
	// the local Delaunay condition needs to be enforced. It is only enforced
	// in one sector, however, that being the angular range defined by
	// fixuptri.
	//
	// This routine also needs to make decisions regarding the "stacking" of
	// triangles.  If the position of
	// the new vertex (the origin of fixuptri) indicates that the vertex before
	// it on the polygon is a reflex vertex, then "stack" the triangle by
	// doing nothing.  (fixuptri is an inverted triangle, which is how stacked
	// triangles are identified.)
	//
	// Otherwise, check whether the vertex before that was a reflex vertex.
	// If so, perform an edge flip, thereby eliminating an inverted triangle
	// (popping it off the stack). The edge flip may result in the creation
	// of a new inverted triangle, depending on whether or not the new vertex
	// is visible to the vertex three edges behind on the polygon.
	//
	// If neither of the two vertices behind the new vertex are reflex
	// vertices, fixuptri and fartri, the triangle opposite it, are not
	// inverted; hence, ensure that the edge between them is locally Delaunay.
	delaunayFixup(fixupTri, leftSide) {
		let nearTri = new Otri();
		let farTri = new Otri();
		let farEdge = new Osub();
		let nearVertex = null;
		let leftVertex = null;
		let rightVertex = null;
		let farVertex = null;

		nearTri = fixupTri.lNext(nearTri);
		farTri = nearTri.sym(farTri);
		// Check if the edge opposite the origin of fixuptri can be flipped.
		if (farTri.Triangle === Mesh.dummyTri()) {
			return;
		}
		farEdge = nearTri.segPivot();
		if (farEdge.Segment !== Mesh.dummySub()) {
			return;
		}
		// Find all the relevant vertices.
		nearVertex = nearTri.apex();
		leftVertex = nearTri.org();
		rightVertex = nearTri.dest();
		farVertex = farTri.apex();

		// Check whether the previous polygon vertex is a reflex vertex.
		if (leftSide === true) {
			if (Primitives.ccw(nearVertex, leftVertex, farVertex) <= 0.0) {
				// leftvertex is a reflex vertex too. Nothing can
				// be done until a convex section is found.
				return;
			}
		} else {
			if (Primitives.ccw(farVertex, rightVertex, nearVertex) <= 0.0) {
				// rightvertex is a reflex vertex too.  Nothing can
				// be done until a convex section is found.
				return;
			}
		}
		if (Primitives.ccw(rightVertex, leftVertex, farVertex) > 0.0) {
			// fartri is not an inverted triangle, and farvertex is not a reflex
			// vertex.  As there are no reflex vertices, fixuptri isn't an
			// inverted triangle, either.  Hence, test the edge between the
			// triangles to ensure it is locally Delaunay.
			if (Primitives.inCircle(leftVertex, farVertex, rightVertex, nearVertex) <= 0.0) {
				return;
			}
			// Not locally Delaunay; go on to an edge flip.
		}
		// else fartri is inverted; remove it from the stack by flipping.
		this.flip(nearTri);
		fixupTri.lPrevSelf(); // Restore the origin of fixuptri after the flip.
		// Recursively process the two triangles that result from the flip.
		this.delaunayFixup(fixupTri, leftSide);
		this.delaunayFixup(farTri, leftSide);
	}

	// Find the intersection of an existing segment and a segment that is being 
	// inserted. Insert a vertex at the intersection, splitting an existing subsegment.
	// The segment being inserted connects the apex of splittri to endpoint2.
	// splitsubseg is the subsegment being split, and must adjoin splittri.
	// Endpoints of the subsegment being split are the origin and
	// destination of splittri.
	//
	// On completion, splittri is a handle having the newly inserted
	// intersection point as its origin, and endpoint1 as its destination.
	segmentIntersection(splitTri, splitSubSeg, endPoint2) {
		let oppSubSeg = new Osub();
		let endPoint1 = null;
		let tOrg = null;
		let tDest = null;
		let leftVertex = null;
		let rightVertex = null;
		let newVertex = null;
		let success = null;

		let ex = 0;
		let ey = 0;
		let tx = 0;
		let ty = 0;
		let etx = 0;
		let ety = 0;
		let split = 0;
		let denom = 0;

		// Find the other three segment endpoints.
		endPoint1 = splitTri.apex();
		tOrg = splitTri.org();
		tDest = splitTri.dest();

        tx = tDest.X - tOrg.X;
		ty = tDest.Y - tOrg.Y;
		ex = endPoint2.X - endPoint1.X;
		ey = endPoint2.Y - endPoint1.Y;
		etx = tOrg.X - -endPoint2.X;
		ety = tOrg.Y - -endPoint2.Y;
		denom = ty * ex - tx * ey;
		if (denom === 0.0) {
			alert('Attempt to find intersection of parallel segments.');
		}
		split = (ey * etx - ex * ety) / denom;
		// Create the new vertex.
		newVertex = new Vertex(0,
			tOrg.X + split * (tDest.X - tOrg.X),
			tOrg.Y + split * (tDest.Y - tOrg.Y),
			this.NExtras,
			splitSubSeg.Segment.BoundaryMarker);
		newVertex.Hash = this.Hash_Vtx++;
		newVertex.Id = newVertex.Hash;
		// Interpolate its attributes.
		for (i = 0; i < this.NExtras; i++) {
			newVertex.Attributes[i] = tOrg.Attributes[i] + split * (tDest.Attributes[i] - tOrg.Attributes[i]);
		}
		this.Vertices.set(newVertex.Hash, newVertex);
		// Insert the intersection vertex.  This should always succeed.
		success = this.insertVertex(newVertex, splitTri, splitSubSeg, false, false);
		if (success !== InsertVertexResult.Successful) {
			alert('Failure to split a segment.');
		}
		// Record a triangle whose origin is the new vertex.
		newVertex.Tri = splitTri;
		if (this.SteinerLeft > 0) {
			this.SteinerLeft--;
		}
		// Divide the segment into two, and correct the segment endpoints.
		splitSubSeg.symSelf();
		oppSubSeg = splitSubSeg.pivot();
		splitSubSeg.dissolve();
		oppSubSeg.dissolve();
		do {
			splitSubSeg.setSegOrg(newVertex);
			splitSubSeg.nextSelf();
		} while (splitSubSeg.Segment !== Mesh.dummySub());
		do {
			oppSubSeg.setSegOrg(newVertex);
			oppSubSeg.nextSelf();
		} while (oppSubSeg.Segment !== Mesh.dummySub());
		// Inserting the vertex may have caused edge flips.  We wish to rediscover
		// the edge connecting endpoint1 to the new intersection vertex.
		this.findDirection(splitTri, endPoint1);
		rightVertex = splitTri.dest();
		leftVertex = splitTri.apex();
		if (leftVertex.X === endPoint1.X && leftVertex.Y === endPoint1.Y) {
			splitTri.oNextSelf();
		} else if (rightVertex.X === endPoint1.X && rightVertex.Y === endPoint1.Y) {
			alert('Topological inconsistency after splitting a segment.');
		}
	}

	// Find the first triangle on the path from one point to another.
	// The return value notes whether the destination or apex of the found
	// triangle is collinear with the two points in question.
	// Finds the triangle that intersects a line segment drawn from the
	// origin of 'searchtri' to the point 'searchpoint', and returns the result
	// in 'searchtri'. The origin of 'searchtri' does not change, even though
	// the triangle returned may differ from the one passed in. This routine
	// is used to find the direction to move in to get from one point to
	// another.
	findDirection(searchTri, searchPoint) {
		let checkTri = new Otri();
		let startVertex = null;
		let leftVertex = null;
		let rightVertex = null;
		let leftCCW = null;
		let rightCCW = null;
		let leftFlag = false;
		let rightFlag = false;

		startVertex = searchTri.org();
		rightVertex = searchTri.dest();
		leftVertex = searchTri.apex();

		// Is 'searchpoint' to the left?
		leftCCW = Primitives.ccw(searchPoint, startVertex, leftVertex);
		leftFlag = leftCCW > 0.0;

		// Is 'searchpoint' to the right?
		rightCCW = Primitives.ccw(startVertex, searchPoint, rightVertex);
		rightFlag = rightCCW > 0.0;
		if (leftFlag === true && rightFlag === true) {
			// 'searchtri' faces directly away from 'searchpoint'. We could go left
			// or right. Ask whether it's a triangle or a boundary on the left.
			checkTri = searchTri.oNext(checkTri);
			if (checkTri.Triangle === Mesh.dummyTri()) {
				leftFlag = false;
			} else {
				rightFlag = false;
			}
		}

		while (leftFlag === true) {
			// Turn left until satisfied.
			searchTri.oNextSelf();
			if (searchTri.Triangle === Mesh.dummyTri()) {
				alert("Unable to find a triangle on path");
			}
			leftVertex = searchTri.apex();
			rightCCW = leftCCW;
			leftCCW = Primitives.ccw(searchPoint, startVertex, leftVertex);
			leftFlag = leftCCW > 0.0;
		}

		while (rightFlag === true) {
			// Turn right until satisfied.
			searchTri.oPrevSelf();
			if (searchTri.Triangle === Mesh.dummyTri()) {
				alert("Unable to find a triangle on path");
			}
			rightVertex = searchTri.dest();
			leftCCW = rightCCW;
			rightCCW = Primitives.ccw(startVertex, searchPoint, rightVertex);
			rightFlag = rightCCW > 0.0;
		}
		if (leftCCW === 0.0) {
			return FindDirectionResult.Leftcollinear;
		} else if (rightCCW === 0.0) {
			return FindDirectionResult.Rightcollinear;
		} else {
			return FindDirectionResult.Within;
		}
	}

	// Construct a mapping from vertices to triangles to improve the speed of 
	// point location for segment insertion.
	// Traverses all the triangles, and provides each corner of each triangle
	// with a pointer to that triangle. Of course, pointers will be overwritten
	// by other pointers because (almost) each vertex is a corner of several
	// triangles, but in the end every vertex will point to some triangle
	// that contains it.
	makeVertexMap() {
		let triOrg = null;
		for (var t of this.Triangles.values()) {
			let tri = new Otri();
			tri.Triangle = t;
			// Check all three vertices of the triangle.
			for (tri.Orient = 0; tri.Orient < 3; tri.Orient++) {
				triOrg = tri.org();
				triOrg.Tri = tri.clone();
			}
		}
	}



	dummyInit() {
		// Initialize the three adjoining triangles to be "outer space." These
		// will be changed by various bonding operations, but their
		// values don't really matter, as long as they can legally be
		// dereferenced.

		let triangle = new Triangle();
		Mesh.setDummyTri(triangle);
		triangle.Neighbors[0].Triangle = triangle;
		triangle.Neighbors[1].Triangle = triangle;
		triangle.Neighbors[2].Triangle = triangle;
		// Set up 'dummysub', the omnipresent subsegment pointed to by any
		// triangle side or subsegment end that isn't attached to a real
		// subsegment.
		let dummySub = new Segment();
		Mesh.setDummySub(dummySub);
		// Initialize the two adjoining subsegments to be the omnipresent
		// subsegment. These will eventually be changed by various bonding
		// operations, but their values don't really matter, as long as they
		// can legally be dereferenced.
		dummySub.SubSegs[0].Segment = dummySub;
		dummySub.SubSegs[1].Segment = dummySub;
		// Initialize the three adjoining subsegments of 'dummytri' to be
		// the omnipresent subsegment.
		triangle.SubSegs[0].Segment = dummySub;
		triangle.SubSegs[1].Segment = dummySub;
		triangle.SubSegs[2].Segment = dummySub;
	}

	// Read the vertices from memory.
	transferNodes(data) {
		this.InVertices = data.Vertices.length;
		this.Mesh_Dim = 2;
		if (this.InVertices < 3) {
			alert("Input must have at least three input vertices.");
		}
		if (data.Vertices[0].Attributes === null) {
			this.NExtras = 0;
		} else {
			this.NExtras = data.Vertices[0].Attributes.length;
		}
		for (var vertex of data.Vertices) {
			vertex.Hash = this.Hash_Vtx++;
			vertex.Id = vertex.Hash;
			this.Vertices.set(vertex.Hash, vertex);
		}
		this.Bounds = data.BoundingBox;
	}
	//Create a new triangle with orientation zero.
	makeTriangle(newOtri) {
		let triangle = new Triangle();
		triangle.Id = this.Hash_Tri++;
		newOtri.Triangle = triangle;
		newOtri.Orient = 0;
		this.Triangles.set(triangle.Id, triangle);
	}
	// Insert a vertex into a Delaunay triangulation, performing flips as necessary 
	// to maintain the Delaunay property.
	// newvertex The point to be inserted.
	// searchtri The triangle to start the search.
	// splitseg Segment to split.
	// segmentflaws Check for creation of encroached subsegments.
	// triflaws Check for creation of bad quality triangles.
	// If a duplicate vertex or violated segment does not prevent the 
	// vertex from being inserted, the return value will be ENCROACHINGVERTEX if 
	// the vertex encroaches upon a subsegment (and checking is enabled), or
	// SUCCESSFULVERTEX otherwise. In either case, 'searchtri' is set to a handle
	// whose origin is the newly inserted vertex.
	// The point 'newvertex' is located. If 'searchtri.triangle' is not NULL,
	// the search for the containing triangle begins from 'searchtri'.  If
	// 'searchtri.triangle' is NULL, a full point location procedure is called.
	// If 'insertvertex' is found inside a triangle, the triangle is split into
	// three; if 'insertvertex' lies on an edge, the edge is split in two,
	// thereby splitting the two adjacent triangles into four. Edge flips are
	// used to restore the Delaunay property. If 'insertvertex' lies on an
	// existing vertex, no action is taken, and the value DUPLICATEVERTEX is
	// returned. On return, 'searchtri' is set to a handle whose origin is the
	// existing vertex.
	//
	// InsertVertex() does not use flip() for reasons of speed; some
	// information can be reused from edge flip to edge flip, like the
	// locations of subsegments.
	// 
	// Param 'splitseg': Normally, the parameter 'splitseg' is set to NULL, 
	// implying that no subsegment should be split. In this case, if 'insertvertex' 
	// is found to lie on a segment, no action is taken, and the value VIOLATINGVERTEX 
	// is returned. On return, 'searchtri' is set to a handle whose primary edge is the 
	// violated subsegment.
	// If the calling routine wishes to split a subsegment by inserting a vertex in it, 
	// the parameter 'splitseg' should be that subsegment. In this case, 'searchtri' 
	// must be the triangle handle reached by pivoting from that subsegment; no point 
	// location is done.
	// 
	// segmentflaws: Flags that indicate whether or not there should
	// be checks for the creation of encroached subsegments. If a newly inserted 
	// vertex encroaches upon subsegments, these subsegments are added to the list 
	// of subsegments to be split if 'segmentflaws' is set.
	// 
	// Param 'triflaws': Flags that indicate whether or not there should be
	// checks for the creation of bad quality triangles. If bad triangles are 
	// created, these are added to the queue if 'triflaws' is set.
	insertVertex(newVertex, searchTri, splitSeg, segmentFlaws, triFlaws) {
		let horiz = new Otri();
		let top = new Otri();
		let botLeft = new Otri();
		let botRight = new Otri();
		let topLeft = new Otri();
		let topRight = new Otri();
		let newBotLeft = new Otri();
		let newBotRight = new Otri();
		let newTopRight = new Otri();
		let botLCasing = new Otri();
		let botRCasing = new Otri();
		let topLCasing = new Otri();
		let topRCasing = new Otri();
		let testTri = new Otri();
		let botLSubSeg = new Osub();
		let botRSubSeg = new Osub();
		let topLSubSeg = new Osub();
		let topRSubSeg = new Osub();
		let brokenSubSeg = new Osub();
		let checkSubSeg = new Osub();
		let rightSubSeg = new Osub();
		let newSubSeg = new Osub();
		let encroached;
		let first = new Vertex();
		let leftVertex = new Vertex();
		let rightVertex = new Vertex();
		let botVertex = new Vertex();
		let topVertex = new Vertex();
		let farVertex = new Vertex();
		let segmentOrg = new Vertex();
		let segmentDest = new Vertex();
		let region = 0;
		let area = 0;
		let success;
		let intersect;
		let doFlip = false;
		let mirrorFlag = false;
		let enq = false;

		if (splitSeg.Segment === null) {
			// Find the location of the vertex to be inserted.  Check if a good
			// starting triangle has already been provided by the caller.
			if (searchTri.Triangle === Mesh.dummyTri()) {
				// Find a boundary triangle.
				horiz.Triangle = Mesh.dummyTri();
				horiz.Orient = 0;
				horiz.symSelf();
				// Search for a triangle containing 'newVertex'.
				intersect = this.Locator.locate(newVertex, horiz);
			} else {
				// Start searching from the triangle provided by the caller.
				searchTri.copy(horiz);
				intersect = this.Locator.preciseLocate(newVertex, horiz, true);
			}
		} else {
			// The calling routine provides the subsegment in which
			// the vertex is inserted.
			searchTri.copy(horiz);
			intersect = LocateResult.OnEdge;
		}

		if (intersect === LocateResult.OnVertex) {
			// There's already a vertex there.  Return in 'searchtri' a triangle
			// whose origin is the existing vertex.
			horiz.copy(searchTri);
			this.Locator.update(horiz);
			return InsertVertexResult.Duplicate;
		}

		if (intersect === LocateResult.OnEdge || intersect === LocateResult.Outside) {
			// The vertex falls on an edge or boundary.
			if (this.CheckSegments === true && splitSeg.Segment === null) {
				// Check whether the vertex falls on a subsegment.
				brokenSubSeg = horiz.segPivot();
				if (brokenSubSeg.Segment !== Mesh.dummySub()) {
					// Return a handle whose primary edge contains the vertex,
					// which has not been inserted.
					horiz.copy(searchTri);
					this.Locator.update(horiz);
					return InsertVertexResult.Violating;
				}
			}
			// Insert the vertex on an edge, dividing one triangle into two (if
			// the edge lies on a boundary) or two triangles into four.
			botRight = horiz.lPrev(botRight);
			botRCasing = botRight.sym(botRCasing);
			topRight = horiz.sym(topRight);
			// Is there a second triangle?  (Or does this edge lie on a boundary?)
			mirrorFlag = topRight.Triangle !== Mesh.dummyTri();
			if (mirrorFlag === true) {
				topRight.lNextSelf();
				topRCasing = topRight.sym(topRCasing);
				this.makeTriangle(newTopRight);
			} else {
				// Splitting a boundary edge increases the number of boundary edges.
				this.HullSize++;
			}
			this.makeTriangle(newBotRight);
			// Set the vertices of changed and new triangles.
			rightVertex = horiz.org();
			leftVertex = horiz.dest();
			botVertex = horiz.apex();
			newBotRight.setOrg(botVertex);
			newBotRight.setDest(rightVertex);
			newBotRight.setApex(newVertex);
			horiz.setOrg(newVertex);
			// Set the region of a new triangle.
			newBotRight.Triangle.Region = botRight.Triangle.Region;
			if (mirrorFlag === true) {
				topVertex = topRight.dest();
				newTopRight.setOrg(rightVertex);
				newTopRight.setDest(topVertex);
				newTopRight.setApex(newVertex);
				topRight.setOrg(newVertex);
				// Set the region of another new triangle.
				newTopRight.Triangle.Region = topRight.Triangle.Region;
			}

			// There may be subsegments that need to be bonded
			// to the new triangle(s).
			if (this.CheckSegments === true) {
				botRSubSeg = botRight.segPivot();
				if (botRSubSeg.Segment !== Mesh.dummySub()) {
					botRight.segDissolve();
					botRSubSeg = newBotRight.segBond(botRSubSeg);
				}
				if (mirrorFlag === true) {
					topRSubSeg = topRight.segPivot();
					if (topRSubSeg.Segment !== Mesh.dummySub()) {
						topRight.segDissolve();
						topRSubSeg = newTopRight.segBond(topRSubSeg);
					}
				}
			}

			// Bond the new triangle(s) to the surrounding triangles.
			botRCasing = newBotRight.bond(botRCasing);
			newBotRight.lPrevSelf();
			botRight = newBotRight.bond(botRight);
			newBotRight.lPrevSelf();
			if (mirrorFlag === true) {
				topRCasing = newTopRight.bond(topRCasing);
				newTopRight.lNextSelf();
				topRight = newTopRight.bond(topRight);
				newTopRight.lNextSelf();
				newBotRight = newTopRight.bond(newBotRight);
			}
			if (splitSeg.Segment !== null) {
				// Split the subsegment into two.
				splitSeg.setDest(newVertex);
				segmentOrg = splitSeg.segOrg();
				segmentDest = splitSeg.segDest();
				splitSeg.symSelf();
				rightSubSeg = splitSeg.pivot();
				this.insertSubSeg(newBotRight, splitSeg.Segment.BoundaryMarker);
				newSubSeg = newBotRight.segPivot();
				newSubSeg.setSegOrg(segmentOrg);
				newSubSeg.setSegDest(segmentDest);
				newSubSeg = splitSeg.bond(newSubSeg);
				newSubSeg.symSelf();
				rightSubSeg = newSubSeg.bond(rightSubSeg);
				splitSeg.symSelf();

				// Transfer the subsegment's boundary marker to the vertex if required.
				if (newVertex.Mark === 0) {
					newVertex.Mark = splitSeg.Segment.BoundaryMarker;
				}
			}
			// Position 'horiz' on the first edge to check for
			// the Delaunay property.
			horiz.lNextSelf();

		} else {
			// Insert the vertex in a triangle, splitting it into three.
			botLeft = horiz.lNext(botLeft);
			botRight = horiz.lPrev(botRight);
			botLCasing = botLeft.sym(botLCasing);
			botRCasing = botRight.sym(botRCasing);
			this.makeTriangle(newBotLeft);
			this.makeTriangle(newBotRight);
			// Set the vertices of changed and new triangles.
			rightVertex = horiz.org();
			leftVertex = horiz.dest();
			botVertex = horiz.apex();
			newBotLeft.setOrg(leftVertex);
			newBotLeft.setDest(botVertex);
			newBotLeft.setApex(newVertex);
			newBotRight.setOrg(botVertex);
			newBotRight.setDest(rightVertex);
			newBotRight.setApex(newVertex);
			horiz.setApex(newVertex);
			// Set the region of the new triangles.
			newBotLeft.Triangle.Region = horiz.Triangle.Region;
			newBotRight.Triangle.Region = horiz.Triangle.Region;
			// There may be subsegments that need to be bonded
			// to the new triangles.
			if (this.Checksegments === true) {
				botLSubSeg = botLeft.segPivot();
				if (botLSubSeg.Segment !== Mesh.dummySub()) {
					botLeft.segDissolve();
					botLSubSeg = newBotLeft.segBond(botLSubSeg);
				}
				botRSubSeg = botRight.segPivot();
				if (botRSubSeg.Segment !== Mesh.dummySub()) {
					botRight.segDissolve();
					botRSubSeg = newBotRight.segBond(botRSubSeg);
				}
			}

			// Bond the new triangles to the surrounding triangles.
			botLCasing = newBotLeft.bond(botLCasing);
			botRCasing = newBotRight.bond(botRCasing);
			newBotLeft.lNextSelf();
			newBotRight.lPrevSelf();
			newBotRight = newBotLeft.bond(newBotRight);
			newBotLeft.lNextSelf();
			newBotLeft = botLeft.bond(newBotLeft);
			newBotRight.lPrevSelf();
			newBotRight = botRight.bond(newBotRight);
		}
		// The insertion is successful by default, unless an encroached
		// subsegment is found.
		success = InsertVertexResult.Successful;
		// Circle around the newly inserted vertex, checking each edge opposite it 
		// for the Delaunay property. Non-Delaunay edges are flipped. 'horiz' is 
		// always the edge being checked. 'first' marks where to stop circling.
		first = horiz.org();
		rightVertex = first;
		leftVertex = horiz.dest();

		// Circle until finished.
		let count1 = 0;
		while (true) {
			count1++;
			// By default, the edge will be flipped.
			doFlip = true;
			if (this.CheckSegments === true) {
				// Check for a subsegment, which cannot be flipped.
				checkSubSeg = horiz.segPivot();
				if (checkSubSeg.Segment !== this.dummySub()) {
					// The edge is a subsegment and cannot be flipped.
					doFlip = false;
				}
			}
			if (doFlip === true) {
				// Check if the edge is a boundary edge.
				top = horiz.sym(top);
				if (top.Triangle === Mesh.dummyTri()) {
					// The edge is a boundary edge and cannot be flipped.
					doFlip = false;
				} else {
					// Find the vertex on the other side of the edge.
					farVertex = top.apex();
					// In the incremental Delaunay triangulation algorithm, any of
					// 'leftvertex', 'rightvertex', and 'farvertex' could be vertices
					// of the triangular bounding box. These vertices must be
					// treated as if they are infinitely distant, even though their
					// "coordinates" are not.
					if (leftVertex.equals(this.InfVertex1) === true ||
						leftVertex.equals(this.InfVertex2) === true ||
						leftVertex.equals(this.InfVertex3) === true) {
						// 'leftvertex' is infinitely distant. Check the convexity of
						// the boundary of the triangulation. 'farvertex' might be
						// infinite as well, but trust me, this same condition should
						// be applied.
						doFlip = Primitives.ccw(newVertex, rightVertex, farVertex) > 0.0;
					} else if (rightVertex.equals(this.InfVertex1) === true ||
						rightVertex.equals(this.InfVertex2) === true ||
						rightVertex.equals(this.InfVertex3) === true) {
						// 'rightvertex' is infinitely distant. Check the convexity of
						// the boundary of the triangulation. 'farvertex' might be
						// infinite as well, but trust me, this same condition should
						// be applied.
						doFlip = Primitives.ccw(farVertex, leftVertex, newVertex) > 0.0;
					} else if (farVertex.equals(this.InfVertex1) === true ||
						farVertex.equals(this.InfVertex2) === true ||
						farVertex.equals(this.InfVertex3) === true) {
						// 'farvertex' is infinitely distant and cannot be inside
						// the circumcircle of the triangle 'horiz'.
						doFlip = false;
					} else {
						// Test whether the edge is locally Delaunay.
						doFlip = Primitives.inCircle(leftVertex, newVertex, rightVertex, farVertex) > 0.0;
					}
					if (doFlip) {
						// Flip the edge 'horiz' by rotating its containing
						// quadrilateral (the two triangles adjacent to 'horiz').
						// Identify the casing of the quadrilateral.
						topLeft = top.lPrev(topLeft);
						topLCasing = topLeft.sym(topLCasing);
						topRight = top.lNext(topRight);
						topRCasing = topRight.sym(topRCasing);
						botLeft = horiz.lNext(botLeft);
						botLCasing = botLeft.sym(botLCasing);
						botRight = horiz.lPrev(botRight);
						botRCasing = botRight.sym(botRCasing);
						// Rotate the quadrilateral one-quarter turn counterclockwise.
						botLCasing = topLeft.bond(botLCasing);
						botRCasing = botLeft.bond(botRCasing);
						topRCasing = botRight.bond(topRCasing);
						topLCasing = topRight.bond(topLCasing);
						if (this.CheckSegments === true) {
							// Check for subsegments and rebond them to the quadrilateral.
							topLSubSeg = topLeft.segPivot();
							botLSubSeg = botLeft.segPivot();
							botRSubseg = botRight.segPivot();
							topRSubSeg = topRight.segPivot();
							if (topLSubSeg.Segment === this.dummySub()) {
								topRight.segDissolve();
							} else {
								topLSubSeg = topRight.segBond(topLSubSeg);
							}
							if (botLSubSeg.Segment === this.dummySub()) {
								topLeft.segDissolve();
							} else {
								botLSubSeg = topLeft.segBond(botLSubSeg);
							}
							if (botRSubSeg.Segment === this.dummySub()) {
								botLeft.segDissolve();
							} else {
								botRSubSeg = botLeft.segBond(botRSubSeg);
							}
							if (topRSubSeg.Segment === this.dummySub()) {
								botRight.segDissolve();
							} else {
								topRSubSeg = botRight.segBond(topRSubSeg);
							}
						}
						// New vertex assignments for the rotated quadrilateral.
						horiz.setOrg(farVertex);
						horiz.setDest(newVertex);
						horiz.setApex(rightVertex);
						top.setOrg(newVertex);
						top.setDest(farVertex);
						top.setApex(leftVertex);

						// Assign region.
						region = Math.min(top.Triangle.Region, horiz.Triangle.Region);
						top.Triangle.Region = region;
						horiz.Triangle.Region = region;

						// On the next iterations, consider the two edges that were exposed (this
						// is, are now visible to the newly inserted vertex) by the edge flip.
						horiz.lPrevSelf();
						leftVertex = farVertex;

					}
				}
			}
			if (doFlip === false) {
				// Look for the next edge around the newly inserted vertex.
				horiz.lNextSelf();
				testTri = horiz.sym(testTri);
				// Check for finishing a complete revolution about the new vertex, or
				// falling outside of the triangulation. The latter will happen when
				// a vertex is inserted at a boundary.
				if (leftVertex.equals(first) === true ||
					testTri.Triangle === Mesh.dummyTri()) {
					// Return a triangle whose origin is the new vertex.
					searchTri = horiz.lNext(searchTri);
					let recentTri = new Otri();
					recentTri = horiz.lNext(recentTri);
					this.Locator.update(recentTri);
					return success;
				}
				// Finish finding the next edge around the newly inserted vertex.
				horiz = testTri.lNext(horiz);
				rightVertex = leftVertex;
				leftVertex = horiz.dest();
			}
		}

	}

	// Create a new subsegment and inserts it between two triangles. Its 
	// vertices are properly initialized.
	// tri: The new subsegment is inserted at the edge 
	// described by this handle
	// subSegMark: The marker 'subsegmark' is applied to the 
	// subsegment and, if appropriate, its vertices
	insertSubSeg(tri, subSegMark) {
		let oppoTri = new Otri();
		let newSubSeg = new Osub();
		let triOrg = new Vertex();
		let triDest = new Vertex();

		triOrg = tri.org();
		triDest = tri.dest();
		// Mark vertices if possible.
		if (triOrg.Mark === 0) {
			triOrg.Mark = subSegMark;
		}
		if (triDest.Mark === 0) {
			triDest.Mark = subSegMark;
		}
		// Check if there's already a subsegment here.
		newSubSeg = tri.segPivot();
		if (newSubSeg.Segment === Mesh.dummySub()) {
			// Make new subsegment and initialize its vertices.
			this.makeSegment(newSubSeg);
			newSubSeg.setOrg(triDest);
			newSubSeg.setDest(triOrg);
			newSubSeg.setSegOrg(triDest);
			newSubSeg.setSegDest(triOrg);
			// Bond new subsegment to the two triangles it is sandwiched between.
			// The facing triangle 'oppotri' might be equal to 'dummytri'
			// (outer space), but the new subsegment is bonded to it all the same.
			newSubSeg = tri.segBond(newSubSeg);
			oppoTri = tri.sym(oppoTri);
			newSubSeg.symSelf();
			newSubSeg = oppoTri.segBond(newSubSeg);
			newSubSeg.Segment.BoundaryMarker = subSegMark;
		} else {
			if (newSubSeg.Segment.BoundaryMarker === 0) {
				newSubSeg.Segment.BoundaryMarker = subSegMark;
			}
		}



	}
	// Create a new subsegment with orientation zero.
	//newsubseg: Reference to the new subseg
	makeSegment(newSubSeg) {
		let seg = new Segment();
		seg.Hash = this.Hash_Seg++;
		newSubSeg.Segment = seg;
		newSubSeg.Orient = 0;

		this.SubSegs.set(seg.Hash, seg);
	}

	// Deallocate space for a triangle, marking it dead.
	triangleDealloc(dyingTriangle) {
		// Mark the triangle as dead. This makes it possible to detect dead 
		// triangles when traversing the list of all triangles.
		Otri.kill(dyingTriangle);
		this.Triangles.delete(dyingTriangle.Id);
	}

	renumber() {
		this.renumberMesh(NodeNumbering.Linear);
	}

	renumberMesh(num) {
		// Don't need to do anything if the nodes are already numbered.
		if (num === this.Numbering) {
			return;
		}
		let id;
		if (num === NodeNumbering.Linear) {
			id = 0;
			for (var node of this.Vertices.values()) {
				node.Id = id++;
			}
		}
		// Remember the current numbering.
		this.Numbering = num;
		// Triangles will always be numbered from 0 to n-1
		id = 0;
		for (var item of this.Triangles.values()) {
			item.Id = id++;
		}
	}

	numberOfEdges() {
		return this.Edges;
	}

	isPolygon() {
		return (this.InSegments > 0);
    }
    
	// Deallocate space for a subsegment, marking it dead.
	subSegDealloc(dyingSubSeg) {
		// Mark the subsegment as dead. This makes it possible to detect dead 
		// subsegments when traversing the list of all subsegments.  
		Osub.kill(dyingSubSeg);
		this.SubSegs.delete(dyingSubSeg.Hash);

	}
}

class ReturnValue {
	constructor(node, lineNumber) {
		this.Node = node;
		this.LineNumber = lineNumber;
	}
}

class BoundingBox {
	constructor(xMin = null, xMax = null, yMin = null, yMax = null) {
		if (xMin !== null) {
			this.XMin = xMin;
			this.YMin = yMin;
			this.XMax = xMax;
			this.YMax = yMax;
		} else {
			this.init();
		}
	}

	init() {
		this.XMin = Number.MAX_VALUE;
		this.YMin = Number.MAX_VALUE;
		this.XMax = -Number.MAX_VALUE;
		this.YMax = -Number.MAX_VALUE;
	}

	update(x, y) {
		this.XMin = Math.min(x, this.XMin);
		this.YMin = Math.min(y, this.YMin);
		this.XMax = Math.max(x, this.XMax);
		this.YMax = Math.max(y, this.YMax);
	}

	recalculate(vertices) {
		this.init();
		for (var i = 0; i < vertices.length; i++) {
			this.update(vertices[i].X, vertices[i].Y);
		}
	}
	width() {
		return Math.abs(this.XMax - this.XMin);
	}

	height() {
		return Math.abs(this.YMax - this.YMin);
    }
    
	// Check if given point is inside bounding box.
	// Return true, if bounding box contains given point.
	contains(pt) {
		return (pt.X >= this.XMin && pt.X <= this.XMax && pt.Y >= this.YMin && pt.Y <= this.YMax);
	}

}

class Incremental {
	triangulate(mesh) {
		this.Mesh = mesh;
		let startTri = new Otri();
		// Create a triangular bounding box.
		this.getBoundingBox();
		for (var v of this.Mesh.Vertices.values()) {
			startTri.Triangle = Mesh.dummyTri();
			let tmp = new Osub();
			if (this.Mesh.insertVertex(v, startTri, tmp, false, false) === InsertVertexResult.Duplicate) {
				v.Type = VertexType.UndeadVertex;
				this.Mesh.Undeads++;
			}
		}
		let rb = this.removeBox();
		return rb;
	}

	// Form an "infinite" bounding triangle to insert vertices into.
	// The vertices at "infinity" are assigned finite coordinates, which are
	// used by the point location routines, but (mostly) ignored by the
	// Delaunay edge flip routines.
	getBoundingBox() {
		let infTri = new Otri();
		let box = this.Mesh.Bounds;
		let vertices = [3];
		let xMin = box.XMin;
		let yMin = box.YMin;
		let xMax = box.XMax;
		let yMax = box.YMax;
		let width = box.width();
		let height = box.height();
		let maxWH = Math.max(width, height);
		if (maxWH === 0.0) {
			maxWH = 1.0;
		}

		this.Mesh.InfVertex1 = new Vertex(0, xMin - 50 * maxWH, yMin - 40 * maxWH, null, null);
		this.Mesh.InfVertex2 = new Vertex(0, xMax + 50 * maxWH, yMin - 40 * maxWH, null, null);
		this.Mesh.InfVertex3 = new Vertex(0, 0.5 * (xMin + xMax), yMax + 60 * maxWH, null, null);

		vertices[0] = this.Mesh.InfVertex1;
		vertices[1] = this.Mesh.InfVertex2;
		vertices[2] = this.Mesh.InfVertex3;
		// Create the bounding box.
		this.Mesh.makeTriangle(infTri);

		infTri.setOrg(this.Mesh.InfVertex1);
		infTri.setDest(this.Mesh.InfVertex2);
		infTri.setApex(this.Mesh.InfVertex3);
		// Link dummytri to the bounding box so we can always find an
		// edge to begin searching (point location) from.
		Mesh.dummyTri().Neighbors[0] = infTri.clone();
	}

	removeBox() {
		let deadTriangle = new Otri();
		let searchEdge = new Otri();
		let checkEdge = new Otri();
		let nextEdge = new Otri();
		let finalEdge = new Otri();
		let dissolveEdge = new Otri();
		let markOrg = new Vertex();

		let hullSize = 0;
		// Find a boundary triangle.
		nextEdge.Triangle = Mesh.dummyTri();
		nextEdge.Orient = 0;
		nextEdge.symSelf();

		// Mark a place to stop.
		finalEdge = nextEdge.lPrev(finalEdge);
		nextEdge.lNextSelf();
		nextEdge.symSelf();

		// Find a triangle (on the boundary of the vertex set) that isn't
		// a bounding box triangle.
		searchEdge = nextEdge.lPrev(searchEdge);
		searchEdge.symSelf();

		// Check whether nextedge is another boundary triangle
		// adjacent to the first one.
		checkEdge = nextEdge.lNext(checkEdge);
		checkEdge.symSelf();
		if (checkEdge.Triangle === Mesh.dummyTri()) {
			// Go on to the next triangle.  There are only three boundary
			// triangles, and this next triangle cannot be the third one,
			// so it's safe to stop here.
			searchEdge.lPrevSelf();
			searchEdge.symSelf();
		}
		// Find a new boundary edge to search from, as the current search
		// edge lies on a bounding box triangle and will be deleted.
		Mesh.dummyTri().Neighbors[0] = searchEdge;
		hullSize = -2;
		while (!nextEdge.equal(finalEdge)) {
			hullSize++;
			dissolveEdge = nextEdge.lPrev(dissolveEdge);
			dissolveEdge.symSelf();

			// Disconnect the bounding box triangle from the mesh triangle.
			dissolveEdge.dissolve();
			deadTriangle = nextEdge.lNext(deadTriangle);
			nextEdge = deadTriangle.sym(nextEdge);
			// Get rid of the bounding box triangle.
			this.Mesh.triangleDealloc(deadTriangle.Triangle);
			// Do we need to turn the corner?
			if (nextEdge.Triangle === Mesh.dummyTri()) {
				// Turn the corner.
				dissolveEdge.copy(nextEdge);
			}
		}
		this.Mesh.triangleDealloc(finalEdge.Triangle);
		return hullSize;
	}
}

class RenderData {
	constructor() {
		this.Vertices = [];
		this.Segments = [];
		this.Triangles = [];
		this.MeshEdges = [];
		this.TrianglePartition = [];

		this.NumberOfInputPoints = 0;
		this.NumberOfRegions = 0;
		this.Bounds = null;
	}

	setInputGeometry(data) {
		// Clear unused buffers
		this.Segments = null;
		this.Triangles = null;
		this.MeshEdges = null;

		let n = data.count();
		let i = 0;
		this.NumberOfInputPoints = n;

		// Copy points
		this.Vertices = [2 * n];
		for (var vertex of data.Vertices) {
			this.Vertices[2 * i] = vertex.X;
			this.Vertices[2 * i + 1] = vertex.Y;
			i++;
		}
		// Copy segments
		n = data.Segments.length;
		if (n > 0) {
			this.Segments = [];
			for (var seg of data.Segments) {
				this.Segments.push(seg.P0);
				this.Segments.push(seg.P1);
			}
		}
		this.Bounds = new BoundingBox(data.BoundingBox.XMin, data.BoundingBox.XMax, data.BoundingBox.YMin, data.BoundingBox.YMax);
	}

	setMesh(mesh) {
		// Clear unused buffers
		this.Segments = null;

		let n = mesh.Vertices.length();
		let i = 0;
		this.NumberOfInputPoints = mesh.Number;
		// Linear numbering of mesh
		mesh.renumber();

		// Copy points
		this.Vertices = [2 * n];
		for (var pt of this.Mesh.Vertices.values()) {
			this.Vertices[2 * i] = pt.X
			this.Vertices[2 * i + 1] = pt.Y
		}

		// Copy segments
		n = mesh.Segments.length();
		if (n > 0 && mesh.isPolygon()) {
			let segments = [];
			for (var seg of this.Mesh.Segments.values()) {
				segments.push(seg.p0());
				segments.push(seg.p1());
			}
			this.Segments = Array.from(segments);
		}

		// Copy edges
		let edges = [];
		let e = new EdgeEnumerator(mesh); 
		while (e.moveNext()) {
			edges.push(e.current.p0());
			edges.push(e.current.p1());
		}
		this.MeshEdges = Array.from(edges);
		if (this.NumberOfRegions > 0) {
			this.TrianglePartition = [mesh.Triangles.length];
		}
		i = 0;

		//Copy Triangles
		let triangles = [3 * mesh.Triangles.length];
		for (var tri of this.Mesh.Triangles.values()) {
			triangles.add(tri.p0());
			triangles.add(tri.p1());
			triangles.add(tri.p1());
		}
		this.Triangles = Array.from(triangles);
		this.Bounds = new BoundingBox(mesh.Bounds.XMin, mesh.Bounds.XMax, mesh.Bounds.YMin, mesh.Bounds.YMax);
	}
}

// Carves holes into the triangulation.
class Carver {
	constructor(mesh) {
		this.Mesh = mesh;
		this.Viri = [];
	}

	// Virally infect all of the triangles of the convex hull that are not 
	// protected by subsegments. Where there are subsegments, set boundary 
	// markers as appropriate.
	infectHull() {
		let hullTri = new Otri();
		let nextTri = new Otri();
		let startTri = new Otri();
		let hullSubSeg = new Osub();
		let hOrg = new Vertex();
		let hDest = new Vertex();

		// Find a triangle handle on the hull.
		hullTri.Triangle = Mesh.dummyTri();
		hullTri.Orient = 0;
		hullTri.symSelf();
		// Remember where we started so we know when to stop.
		hullTri.copy(startTri);
		// Go once counterclockwise around the convex hull.
		do {
			// Ignore triangles that are already infected.
			if (!hullTri.isInfected()) {
				// Is the triangle protected by a subsegment?
				hullSubSeg = hullTri.segPivot(hullSubSeg);
				if (hullSubSeg.Segment === Mesh.dummySub()) {
					// The triangle is not protected; infect it.
					if (!hullTri.isInfected()) {
						hullTri.infect();
						this.Viri.push(hullTri.Triangle);
					}
				} else {
					// The triangle is protected; set boundary markers if appropriate.
					if (hullSubSeg.Segment.BoundaryMarker === 0) {
						hullSubSeg.Segment.BoundaryMarker = 1;
						hOrg = hullTri.org();
						hDest = hullTri.dest();
						if (hOrg.Mark === 0) {
							hOrg.Mark = 1;
						}
						if (hDest.Mark === 0) {
							hDest.Mark = 1;
						}
					}

				}
			}
			// To find the next hull edge, go clockwise around the next vertex.
			hullTri.lNextSelf();
			nextTri = hullTri.oPrev(nextTri);
			while (nextTri.Triangle !== Mesh.dummyTri()) {
				nextTri.copy(hullTri);
				nextTri = hullTri.oPrev(nextTri);
			}
		} while (!hullTri.equal(startTri));

    }
    
	// Spread the virus from all infected triangles to any neighbors not 
	// protected by subsegments. Delete all infected triangles.
	// This is the procedure that actually creates holes and concavities.
	//
	// This procedure operates in two phases. The first phase identifies all
	// the triangles that will die, and marks them as infected. They are
	// marked to ensure that each triangle is added to the virus pool only
	// once, so the procedure will terminate.
	//
	// The second phase actually eliminates the infected triangles. It also
	// eliminates orphaned vertices.
	plague() {
		let testTri = new Otri();
		let neighbor = new Otri();
		let neighborSubSeg = new Osub();
		let testVertex = new Vertex();
		let nOrg = new Vertex();
		let nDest = new Vertex();
		let killOrg = false;
		// Loop through all the infected triangles, spreading the virus to
		// their neighbors, then to their neighbors' neighbors.
		for (var i = 0; i < this.Viri.length; i++) {
			testTri.Triangle = this.Viri[i];
			// A triangle is marked as infected by messing with one of its pointers
			// to subsegments, setting it to an illegal value.  Hence, we have to
			// temporarily uninfect this triangle so that we can examine its
			// adjacent subsegments.
			testTri.uninfect();
			// Check each of the triangle's three neighbors.

			for (testTri.Orient = 0; testTri.Orient < 3; testTri.Orient++) {
				// Find the neighbor.
				neighbor = testTri.sym(neighbor);
				// Check for a subsegment between the triangle and its neighbor.
				neighborSubSeg = testTri.segPivot(neighborSubSeg);
				// Check if the neighbor is nonexistent or already infected.
				if (neighbor.Triangle === Mesh.dummyTri() || neighbor.isInfected()) {
					if (neighborSubSeg.Segment !== Mesh.dummySub()) {
						// There is a subsegment separating the triangle from its
						// neighbor, but both triangles are dying, so the subsegment
						// dies too.
						this.Mesh.subSegDealloc(neighborSubSeg.Segment);
						if (neighbor.Triangle !== Mesh.dummytri()) {
							// Make sure the subsegment doesn't get deallocated again
							// later when the infected neighbor is visited.
							neighbor.uninfect();
							neighbor.segDissolve();
							neighbor.infect();
						}
					}
				} else {
					// The neighbor exists and is not infected.
					if (neighborSubSeg.Segment === Mesh.dummySub()) {
						// There is no subsegment protecting the neighbor, so
						// the neighbor becomes infected.
						neighbor.infect();
						// Ensure that the neighbor's neighbors will be infected.
						this.Viri.push(neighbor.Triangle);
					} else {
						// The neighbor is protected by a subsegment.
						// Remove this triangle from the subsegment.
						neighborSubSeg.triDissolve();
						// The subsegment becomes a boundary.  Set markers accordingly.
						if (neighborSubSeg.Segment.BoundaryMarker === 0) {
							neighborSubSeg.Segment.BoundaryMarker = 1;
						}
						nOrg = neighbor.org();
						nDest = neighbor.dest();
						if (nOrg.Mark === 0) {
							nOrg.Mark = 1;
						}
						if (nDest.Mark === 0) {
							nDest.Mark = 1;
						}
					}

				}

			}
			// Remark the triangle as infected, so it doesn't get added to the
			// virus pool again.
			testTri.infect();
		}
		for (var virus of this.Viri) {
			testTri.Triangle = virus;
			// Check each of the three corners of the triangle for elimination.
			// This is done by walking around each vertex, checking if it is
			// still connected to at least one live triangle.
			for (testTri.Orient = 0; testTri.Orient < 3; testTri.Orient++) {
				testVertex = testTri.org();
				// Check if the vertex has already been tested.
				if (testVertex !== null) {
					killOrg = true;
					// Mark the corner of the triangle as having been tested.
					testTri.setOrg(null);
					// Walk counterclockwise about the vertex.
					neighbor = testTri.oNext(neighbor);
					// Stop upon reaching a boundary or the starting triangle.
					while (neighbor.Triangle !== Mesh.dummyTri() && !neighbor.equal(testTri)) {
						if (neighbor.isInfected()) {
							// Mark the corner of this triangle as having been tested.
							neighbor.setOrg(null);
						} else {
							// A live triangle.  The vertex survives.
							killOrg = false;
						}
						// Walk counterclockwise about the vertex.
						neighbor.oNextSelf();
					}
					// If we reached a boundary, we must walk clockwise as well.
					if (neighbor.Triangle === Mesh.dummyTri()) {
						// Walk clockwise about the vertex.
						neighbor = testTri.oPrev(neighbor);
						// Stop upon reaching a boundary.
						while (neighbor.Triangle !== Mesh.dummyTri()) {
							if (neighbor.isInfected()) {
								// Mark the corner of this triangle as having been tested.
								neighbor.setOrg(null);
							} else {
								// A live triangle.  The vertex survives.
								killOrg = false;
							}
							// Walk clockwise about the vertex.
							neighbor.oPrevSelf();
						}
					}
					if (killOrg) {
						// Deleting vertex
						testVertex.Type = VertexType.UndeadVertex;
						this.Mesh.Undeads++;
					}
				}
			}
			// Record changes in the number of boundary edges, and disconnect
			// dead triangles from their neighbors.
			for (testTri.Orient = 0; testTri.Orient < 3; testTri.Orient++) {
				neighbor = testTri.sym(neighbor);
				if (neighbor.Triangle === Mesh.dummyTri()) {
					// There is no neighboring triangle on this edge, so this edge
					// is a boundary edge. This triangle is being deleted, so this
					// boundary edge is deleted.
					this.Mesh.HullSize--;
				} else {
					// Disconnect the triangle from its neighbor.
					neighbor.dissolve();
					// There is a neighboring triangle on this edge, so this edge
					// becomes a boundary edge when this triangle is deleted.
					this.Mesh.HullSize++;
				}
			}
			// Return the dead triangle to the pool of triangles.
			this.Mesh.triangleDealloc(testTri.Triangle);
		}
		// Empty the virus pool.
		this.Viri = [];
	}

	// Find the holes and infect them. Find the area constraints and infect 
	// them. Infect the convex hull. Spread the infection and kill triangles. 
	// Spread the area constraints.
	carveHoles() {
		let searchTri = new Otri();
		let searchOrg = new Vertex();
		let seacrhDest = new Vertex();
		let intersect;
		let regionTris = [];

		if (!this.Mesh.Behavior.Convex) {
			// Mark as infected any unprotected triangles on the boundary.
			// This is one way by which concavities are created.
			this.infectHull();
		}

		// Infect each triangle in which a hole lies.
		for (var hole of this.Mesh.Holes) {
			// Ignore holes that aren't within the bounds of the mesh.
			if (this.Mesh.Bounds.contains(hole)) {
				// Start searching from some triangle on the outer boundary.
				searchTri.Triangle = Mesh.dummyTri();
				searchTri.Orient = 0;
				searchTri.symSelf();
				// Ensure that the hole is to the left of this boundary edge;
				// otherwise, locate() will falsely report that the hole
				// falls within the starting triangle.
				searchOrg = searchTri.org();
				seacrhDest = searchTri.dest();

				if (Primitives.ccw(searchOrg, seacrhDest, hole) > 0.0) {
					// Find a triangle that contains the hole.
					intersect = this.Mesh.Locator.locate(hole, searchTri);
					if (intersect != LocateResult.Outside && !searchTri.isInfected()) {
						// Infect the triangle. This is done by marking the triangle
						// as infected and including the triangle in the virus pool. 
						searchTri.infect();
						this.Viri.push(searchTri.Triangle);
					}
				}


			}
		}
		// Find all the regions before we carve the holes, because locate() won't
		// work when the triangulation is no longer convex. 
		if (this.Mesh.Regions.length > 0) {
			let i = 0;
			regionTris = [this.Mesh.Regions.Lenght];
			// Find the starting triangle for each region.
			for (var region of this.Mesh.Regions) {
				regionTris[i] = Mesh.dummyTri();
				// Ignore region points that aren't within the bounds of the mesh.
				if (this.Mesh.Bounds.contains(region.Point)) {
					// Start searching from some triangle on the outer boundary.
					searchTri.Triangle = Mesh.dummyTri();
					searchTri.Orient = 0;
					searchTri.symSelf();
					// Ensure that the region point is to the left of this boundary
					// edge; otherwise, locate() will falsely report that the
					// region point falls within the starting triangle.
					searchOrg = searchTri.org();
					seacrhDest = searchTri.dest();
					if (Primitives.ccw(searchOrg, seacrhDest, region.Point) > 0.0) {
						// Find a triangle that contains the region point.
						intersect = this.Mesh.Locator.locate(region.Point, searchTri);
						if (intersect != LocateResult.Outside && !searchTri.isInfected()) {
							// Record the triangle for processing after the
							// holes have been carved.
							regionTris[i] = searchTri.Triangle;
							regionTris[i].Region = region.Id;
						}

					}

				}
				i++;
			}
		}

		if (this.Viri.length > 0) {
			// Carve the holes and concavities.
			this.plague();
		}
		if (regionTris != null) {
			let iterator = new RegionIterator(this.Mesh);
			for (var i = 0; i < regionTris.Length; i++) {
				if (regionTris[i] !== Mesh.dummyTri()) {
					// Make sure the triangle under consideration still exists.
					// It may have been eaten by the virus. 
					if (!Otri.isDead(regionTris[i])) {
						iterator.process(regionTris[i]);
					}
				}
			}
		}
		// Free up memory (virus pool should be empty anyway).
		this.Viri = [];
	}

}

class MeshFileReader {
	constructor() {
		this.StartIndex = 0;
	}

	readVertices(lines, lineIndex, numOfVertices, node) {
		let verticesIndex = 0;
		var line = lineIndex;
		for (; line < lines.length; line++) {
			var columns = lines[line].trim().split(/\s+/);
			if (columns[0].charAt(0) === '#') {
				continue;
			}
			verticesIndex++;
			if (verticesIndex > numOfVertices) {
				return --line;
			}
			var id = parseFloat(columns[0]);
			var x = parseFloat(columns[1]);
			var y = parseFloat(columns[2]);
			var attributes = parseInt(columns[3]);
			var boundaryMarker = parseFloat(columns[4]);
			var vertex = new Vertex(id, x, y, attributes, boundaryMarker);
			node.addVertex(vertex);
			if (verticesIndex == 1) {
				this.StartIndex = id;
			}
		}
		return --line;
	}

	getNode(lines, startIndex) {
		let node;
		var readFirstLine = false;
		var line = startIndex;
		for (; line < lines.length; line++) {
			var columns = lines[line].trim().split(/\s+/);
			if (columns[0].charAt(0) === '#') {
				continue;
			}
			if (readFirstLine === false) {
				readFirstLine = true;
				var numOfVertices = parseInt(columns[0]);
				var dimension = parseInt(columns[1]);
				var numOfAttributes = parseInt(columns[2]);
				var numOfBoundaryMakers = parseInt(columns[3]);
				node = new Node(numOfVertices, dimension, numOfAttributes, numOfBoundaryMakers);
			} else {
				line = this.readVertices(lines, line, numOfVertices, node);
				node.BoundingBox.recalculate(node.Vertices);
				break;
			}
		}
		return new ReturnValue(node, line);
	}

	readNodeFile(file, onLoadCallback) {
		let reader = new FileReader();
		var context = this;
		reader.onload = function(progressEvent) {
			var lines = this.result.split('\n');
			let node = context.getNode(lines, 0).Node;
			let inputGeometry = new InputGeometry(node.NumOfVertices,
				node.Dimension,
				node.NumOfAttributes,
				node.NumOfBoundaryMakers,
				node.Vertices,
				node.BoundingBox,
				null,
				null,
				null,
				null,
				null,
				null);
			onLoadCallback(inputGeometry);
		};
		reader.readAsText(file);
	}

    readPolyFile(file, nodeFile, onLoadCallback) {
		var context = this;
		if (nodeFile !== null) {
			let reader = new FileReader();
			reader.onload = function(progressEvent) {
				var nodeLines = this.result.split('\n');
				context.readPolyFile2(file, nodeLines, onLoadCallback);
			};
			reader.readAsText(nodeFile);
		} else {
			context.readPolyFile2(file, null, onLoadCallback);
		}
	}

	readPolyFile2(file, nodeLines, onLoadCallback) {
		let reader = new FileReader();
		var context = this;
		reader.onload = function(progressEvent) {
			var lines = this.result.split('\n');
			let poly = context.getPoly(lines, nodeLines);
			let inputGeometry = new InputGeometry(poly.Node.NumOfVertices,
				poly.Node.Dimension,
				poly.Node.NumOfAttributes,
				poly.Node.NumOfBoundaryMakers,
				poly.Node.Vertices,
				poly.Node.BoundingBox,
				poly.NumOfSegments,
				poly.NumOfBoundaryMakers,
				poly.NumberOfHoles,
				poly.NumberOfRegionalAttrs,
				poly.Segments,
				poly.Holes);
			onLoadCallback(inputGeometry);
		};
		reader.readAsText(file);
	}

	getPoly(lines, nodeLines) {
		let readFirstLine = false;
		let segmentHeader = false;
		let holesHeader = false;
		let node;
		let poly;
		var context = this;
		for (var line = 0; line < lines.length; line++) {
			var columns = lines[line].trim().split(/\s+/);
			if (columns[0].charAt(0) === '#') {
				continue;
			}
			if (readFirstLine === false) {
				readFirstLine = true;
				var numOfVertices = parseInt(columns[0]);
				let retVal;
				if (numOfVertices === 0) {
					if (nodeLines === null) {
						alert('.node file must be supplied');
					}
					retVal = this.getNode(nodeLines, 0);
					node = retVal.Node;
					segmentHeader = true;
				} else {
					retVal = this.getNode(lines, line);
					node = retVal.Node;
					line = retVal.LineNumber;
					segmentHeader = true;
				}
				continue;
			}
			if (segmentHeader === true) {
				var numOfSegments = parseInt(columns[0]);
				var numOfBoundaryMakers = parseInt(columns[1]);
				poly = new Poly(node, numOfSegments, numOfBoundaryMakers);
				segmentHeader = false;
				line++;
				line = this.addSegments(poly, lines, line, numOfSegments);
				holesHeader = true;
				continue;
			}
			if (holesHeader === true) {
				let numberOfHoles = parseInt(columns[0]);
				poly.NumberOfHoles = numberOfHoles;
				holesHeader = false;
				line++;
				line = this.addHoles(poly, lines, line, numberOfHoles);
				continue;
			}

		}
		return poly;
	}

	addSegments(poly, lines, startIndex, segmentCount) {
		let count = segmentCount + startIndex;
		for (var line = startIndex; line < count;) {
			var columns = lines[line].trim().split(/\s+/);
			if (columns[0].charAt(0) === '#') {
				line++;
				count++;
				continue;
			}
			let segmentNum = parseInt(columns[0]);
			let endPoint1 = parseInt(columns[1]) - this.StartIndex;
			let endPoint2 = parseInt(columns[2]) - this.StartIndex;
			let boundaryMarker = parseInt(columns[3]);
			let segment = new Edge(endPoint1, endPoint2, boundaryMarker);
			poly.addSegment(segment);
			line++;
		}
		line--;
		return line;
	}

	addHoles(poly, lines, startIndex, holesCount) {
		let count = holesCount + startIndex;
		for (var line = startIndex; line < count;) {
			var columns = lines[line].trim().split(/\s+/);
			if (columns[0].charAt(0) === '#') {
				line++;
				count++;
				continue;
			}
			let id = parseInt(columns[0]);
			let x = parseFloat(columns[1]);
			let y = parseFloat(columns[2]);
			let hole = new Hole(id, x, y);
			poly.addHole(hole);
			line++;
		}
		line--;
		return line;
	}

}

class EdgeEnumerator {
	constructor() {
		this._data = [1, 2, 3, 4];
	}

	[Symbol.iterator]() {
		var index = -1;
		var data = this._data;

		return {
			next: () => ({
				value: data[++index],
				done: !(index in data)
			})
		};
	};
}
// Iterates the region a given triangle belongs to and applies an action
// to each connected trianlge in that region. Default action is to set the 
// region id.
class RegionIterator {
	constructor(mesh) {
		this.Mesh = mesh;
		this.Viri = [];
	}
	// Process all triangles connected to given triangle and apply given action.

	process(triangle) {
		if (triangle !== Mesh.dummyTri()) {
			// Make sure the triangle under consideration still exists.
			// It may have been eaten by the virus.
			if (!Otri.isDead(triangle)) {
				// Put one triangle in the virus pool.
				triangle.Infected = true;
				this.Viri.push(triangle);
				// Apply one region's attribute and/or area constraint.
				this.processRegion(triangle);
				// The virus pool should be empty now.
			}

		}
		// Free up memory (virus pool should be empty anyway).
		this.Viri = [];
	}

	// Spread regional attributes and/or area constraints (from a .poly file) 
	// throughout the mesh.
	// This procedure operates in two phases. The first phase spreads an
	// attribute and/or an area constraint through a (segment-bounded) region.
	// The triangles are marked to ensure that each triangle is added to the
	// virus pool only once, so the procedure will terminate.
	//
	// The second phase uninfects all infected triangles, returning them to
	// normal.
	processRegion(triangle) {
		let testTri = new Otri();
		let neighbor = new Otri();
		let neighborSubSeg = new Otri();
		let behavior = this.Mesh.Behavior;

		// Loop through all the infected triangles, spreading the attribute
		// and/or area constraint to their neighbors, then to their neighbors'
		// neighbors.
		for (var i = 0; i < this.Viri.length; i++) {
			testTri.Triangle = this.Viri[i];
			// A triangle is marked as infected by messing with one of its pointers
			// to subsegments, setting it to an illegal value.  Hence, we have to
			// temporarily uninfect this triangle so that we can examine its
			// adjacent subsegments.
			testTri.uninfect();
			// Apply function.
			testTri.Triangle.region = triangle.Region;
			// Check each of the triangle's three neighbors.
			for (testTri.Orient = 0; testTri.Orient < 3; testTri.Orient++) {
				// Find the neighbor.
				neighbor = testTri.sym(neighbor);
				// Check for a subsegment between the triangle and its neighbor.
				neighborSubSeg = testTri.segPivot(neighborSubSeg);
				// Make sure the neighbor exists, is not already infected, and
				// isn't protected by a subsegment.
				if (neighbor.Triangle !== Mesh.dummyTri() && !neighbor.isInfected() && neighborSubSeg.Segment == Mesh.dummySub()) {
					// Infect the neighbor.
					neighbor.infect();
					// Ensure that the neighbor's neighbors will be infected.
					this.Viri.push(neighbor.Triangle);
				}
			}
			// Remark the triangle as infected, so it doesn't get added to the
			// virus pool again.
			testTri.infect();
		}
		// Uninfect all triangles.
		for (var virus of this.Viri) {
			virus.Infected = false;
		}
		// Empty the virus pool.
		this.Viri = [];
	}
}

var LocateResult = {
	InTriangle: 0,
	OnEdge: 1,
	OnVertex: 2,
	Outside: 3
};
var InsertVertexResult = {
	Successful: 0,
	Encroaching: 1,
	Violating: 2,
	Duplicate: 3
};

var FindDirectionResult = {
	Within: 0,
	Leftcollinear: 1,
	Rightcollinear: 2
};

var VertexType = {
	InputVertex: 0,
	SegmentVertex: 1,
	FreeVertex: 2,
	DeadVertex: 3,
	UndeadVertex: 4
};

var NodeNumbering = {
	None: 0,
	Linear: 1,
	CuthillMcKee: 2
};
