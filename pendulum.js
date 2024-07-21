"use strict";
// var G = 1.2; // gravitational acceleration
const M = 1.0; // mass
const L = 1.0; // length
const dtMax = 30.0; // ms
const tailMax = 400; // tail length

const barWidth = 0.04;
const barLength = 0.23;
const massRadius = 0.035;
const tailThickness = 0.012;

// To draw a filled rectangle centered at the given coordinates
const quad = new Float32Array([-1, -1, +1, -1, -1, +1, +1, +1]);

const massShader = {
    vert: `
precision mediump float;
attribute vec2 a_point; //position of a vertex
uniform   vec2 u_center; //center position of the mass
uniform   vec2 u_aspect; // aspect ratio of viewport
varying   vec2 v_point; // variable to pass vertex position to fragment shader
void main() {
    v_point = a_point;
    gl_Position = vec4(a_point * ${massRadius} / u_aspect + u_center, 0, 1);
}`,
    frag: `
precision mediump float;
uniform vec2 u_aspect;
uniform vec3 u_color;
varying vec2 v_point; // variable from vertex shader
void main() {
    //calculate distance from the center to the current fragment
    float dist = distance(vec2(0, 0), v_point);
    // smooth transition function to create a gradient effect
    float v = smoothstep(1.0, 0.9, dist);
    gl_FragColor = vec4(u_color, v);
}`,
};

const barShader = {
    vert: `
precision mediump float;
attribute vec2  a_point; //position of a vertex
uniform   float u_angle; // rotation angle of the bar
uniform   vec2  u_attach; // attachment point of the bar
uniform   vec2  u_aspect; //aspect ratio of the viewport
void main() {
    // 2x2 rotation matrix 
    mat2 rotate = mat2(+cos(u_angle), +sin(u_angle),
                       -sin(u_angle), +cos(u_angle));
    // scales and translate the vertices to match the bar's width
    vec2 pos = rotate * (a_point * vec2(1, ${barWidth}) + vec2(1, 0));
    gl_Position = vec4((pos * ${barLength} / u_aspect + u_attach), 0, 1);
}`,
    frag: `
precision mediump float;
uniform vec3 u_color;
void main() {
    gl_FragColor = vec4(u_color, 1);
}`,
};

const tailShader = {
    vert: `
precision mediump float;
attribute vec2  a_point; // position of a vertex
attribute float a_alpha; // alpha value for transparency
uniform   vec2  u_aspect; // aspect ratio of the viewport
varying   float v_alpha; // pass alpha value to fragment shader
void main() {
    v_alpha = a_alpha;
    gl_Position = vec4(a_point * vec2(1, -1) / u_aspect, 0, 1);
}`,
    frag: `
precision mediump float;
uniform vec3  u_color; // color of the tail
uniform float u_cutoff; // cutoff value for transparency
varying float v_alpha; // from vertex shader indicating alpha value
void main() {
    float icutoff = 1.0 - u_cutoff;
    // calculate the trasparency of the color based on the age
    gl_FragColor = vec4(u_color, max(0.0, v_alpha - u_cutoff) / icutoff);
}`,
};

function deriviative(a1, a2, p1, p2) {
    let ml2 = M * L * L;
    let cos12 = Math.cos(a1 - a2);
    let sin12 = Math.sin(a1 - a2);
    //rate of change of angle a1
    let da1 = 6 / ml2 * (2 * p1 - 3 * cos12 * p2) / (16 - 9 * cos12 * cos12);
    //rate of change of angle a2 
    let da2 = 6 / ml2 * (8 * p2 - 3 * cos12 * p1) / (16 - 9 * cos12 * cos12); 
     // rate of change of momentum p1
    let dp1 = ml2 / -2 * (+da1 * da2 * sin12 + 3 * G / L * Math.sin(a1));
    // rate of change of momentum p2
    let dp2 = ml2 / -2 * (-da1 * da2 * sin12 + 3 * G / L * Math.sin(a2)); 
    return [da1, da2, dp1, dp2];
}

// Update pendulum by timestep
function rk4(k1a1, k1a2, k1p1, k1p2, dt) {
    // Compute the derivatives at the initial state (k1a1, k1a2, k1p1, k1p2)
    let [k1da1, k1da2, k1dp1, k1dp2] = deriviative(k1a1, k1a2, k1p1, k1p2);
    // Compute the state at the midpoint using the first stage derivatives
    let k2a1 = k1a1 + k1da1 * dt / 2;
    let k2a2 = k1a2 + k1da2 * dt / 2;
    let k2p1 = k1p1 + k1dp1 * dt / 2;
    let k2p2 = k1p2 + k1dp2 * dt / 2;
    // Compute the derivatives at the midpoint state (k2a1, k2a2, k2p1, k2p2)
    let [k2da1, k2da2, k2dp1, k2dp2] = deriviative(k2a1, k2a2, k2p1, k2p2);
    // Compute the state at the midpoint using the second stage derivatives
    let k3a1 = k1a1 + k2da1 * dt / 2;
    let k3a2 = k1a2 + k2da2 * dt / 2;
    let k3p1 = k1p1 + k2dp1 * dt / 2;
    let k3p2 = k1p2 + k2dp2 * dt / 2;
    // Compute the derivatives at the midpoint state (k3a1, k3a2, k3p1, k3p2)
    let [k3da1, k3da2, k3dp1, k3dp2] = deriviative(k3a1, k3a2, k3p1, k3p2);
    // Compute the state at the end of the interval using the third stage derivatives:
    let k4a1 = k1a1 + k3da1 * dt;
    let k4a2 = k1a2 + k3da2 * dt;
    let k4p1 = k1p1 + k3dp1 * dt;
    let k4p2 = k1p2 + k3dp2 * dt;
    // Combine the derivatives from the four stages to get the final state
    let [k4da1, k4da2, k4dp1, k4dp2] = deriviative(k4a1, k4a2, k4p1, k4p2);
    // Returns an array containing the updated values
    return [
        k1a1 + (k1da1 + 2*k2da1 + 2*k3da1 + k4da1) * dt / 6,
        k1a2 + (k1da2 + 2*k2da2 + 2*k3da2 + k4da2) * dt / 6,
        k1p1 + (k1dp1 + 2*k2dp1 + 2*k3dp1 + k4dp1) * dt / 6,
        k1p2 + (k1dp2 + 2*k2dp2 + 2*k3dp2 + k4dp2) * dt / 6
    ];
}

function history(n) {
    let h = {
        i: 0, //Index to track the current position in the circular buffer
        length: 0, // The current number of elements stored in the history (up to n).
        v: new Float32Array(n * 2), // Float32Array of size n * 2 to store the sine and cosine sums of the angles
        // Adds a new entry to the history
        push: function(a1, a2) { 
            h.v[h.i * 2 + 0] = Math.sin(a1) + Math.sin(a2);
            h.v[h.i * 2 + 1] = Math.cos(a1) + Math.cos(a2);
            h.i = (h.i + 1) % n;
            if (h.length < n)
                h.length++;
        },
        // Applies a function f to each pair of consecutive entries in the history
        visit: function(f) {
            for (let j = h.i + n - 2; j > h.i + n - h.length - 1; j--) {
                let a = (j + 1) % n;
                let b = (j + 0) % n;
                f(h.v[a * 2], h.v[a * 2 + 1], h.v[b * 2], h.v[b * 2 + 1]);
            }
        }
    };
    return h;
}

function normalize(v0, v1) {
    let d = Math.sqrt(v0 * v0 + v1 * v1);
    return [v0 / d, v1 / d];
}

function sub(a0, a1, b0, b1) {
    return [a0 - b0, a1 - b1];
}

function add(a0, a1, b0, b1) {
    return [a0 + b0, a1 + b1];
}

function dot(ax, ay, bx, by) {
    return ax * bx + ay * by;
}

/* Convert tail line into a triangle strip.
 * https://forum.libcinder.org/topic/smooth-thick-lines-using-geometry-shader
 */
function polyline(hist, poly) {
    const w = tailThickness;
    let i = -1;
    let x0, y0;
    let xf, yf;
    hist.visit(function(x1, y1, x2, y2) {
        if (++i === 0) { // For the first point (x1, y1), we initialize the polyline
            let [lx, ly] = sub(x2, y2, x1, y1);
            let [nx, ny] = normalize(-ly, lx);
            poly[0] = x1 + w * nx;
            poly[1] = y1 + w * ny;
            poly[2] = x1 - w * nx;
            poly[3] = y1 - w * ny;
        } else { // For subsequent points, we handle each segment between points (x1, y1) and (x2, y2)
            let [ax, ay] = sub(x1, y1, x0, y0);
            [ax, ay] = normalize(ax, ay);
            let [bx, by] = sub(x2, y2, x1, y1);
            [bx, by] = normalize(bx, by);
            let [tx, ty] = add(ax, ay, bx, by);
            [tx, ty] = normalize(tx, ty);
            let [mx, my] = [-ty, tx];
            let [lx, ly] = sub(x1, y1, x0, y0);
            let [nx, ny] = normalize(-ly, lx);
            let len = Math.min(w, w / dot(mx, my, nx, ny));
            poly[i * 4 + 0] = x1 + mx * len;
            poly[i * 4 + 1] = y1 + my * len;
            poly[i * 4 + 2] = x1 - mx * len;
            poly[i * 4 + 3] = y1 - my * len;
        }
        x0 = x1;
        y0 = y1;
        xf = x2;
        yf = y2;
    });
    // After visiting all points, handle the final point (xf, yf)
    let [lx, ly] = sub(xf, yf, x0, y0);
    let [nx, ny] = normalize(-ly, lx);
    i++;
    poly[i * 4 + 0] = xf + w * nx;
    poly[i * 4 + 1] = yf + w * ny;
    poly[i * 4 + 2] = xf - w * nx;
    poly[i * 4 + 3] = yf - w * ny;
}


// Compiling the given vertex and fragment shader source code into a program.
function compile(gl, vert, frag) {
    let v = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(v, vert);
    let f = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(f, frag);
    gl.compileShader(v);
    if (!gl.getShaderParameter(v, gl.COMPILE_STATUS))
        throw new Error(gl.getShaderInfoLog(v));
    gl.compileShader(f);
    if (!gl.getShaderParameter(f, gl.COMPILE_STATUS))
        throw new Error(gl.getShaderInfoLog(f));
    let p = gl.createProgram();
    gl.attachShader(p, v);
    gl.attachShader(p, f);
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS))
        throw new Error(gl.getProgramInfoLog(p));
    gl.deleteShader(v);
    gl.deleteShader(f);
    let result = {
        program: p
    };
    let nattrib = gl.getProgramParameter(p, gl.ACTIVE_ATTRIBUTES);
    for (let a = 0; a < nattrib; a++) {
        let name = gl.getActiveAttrib(p, a).name;
        result[name] = gl.getAttribLocation(p, name);
    }
    let nuniform = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
    for (let u = 0; u < nuniform; u++) {
        let name = gl.getActiveUniform(p, u).name;
        result[name] = gl.getUniformLocation(p, name);
    }
    return result;
};

// Create a new, random double pendulum
function pendulum({
    tailColor = [0, 0, 1],
    massColor = [0, 0, 0],
    init = null
} = {}) {
    let tail = new history(tailMax); // Create a new history object for the tail
    let a1, a2, p1, p2; // Variables for angles and momenta
    if (init) {
        [a1, a2, p1, p2] = init;
    } else {
        a1 = Math.random() * Math.PI / 2 + Math.PI * 3 / 4;
        a2 = Math.random() * Math.PI / 2 + Math.PI * 3 / 4;
        p1 = 0.0;
        p2 = 0.0;
    }

    return {
        tailColor: tailColor,
        massColor: massColor,
        tail: tail,
        state: function() {
            return [a1, a2, p1, p2]; //angle (a1, a2) and momentum (p1, p2)
        },
        positions: function() {
            let x1 = +Math.sin(a1);
            let y1 = -Math.cos(a1);
            let x2 = +Math.sin(a2) + x1;
            let y2 = -Math.cos(a2) + y1;
            return [x1, y1, x2, y2];
        },
        step: function(dt) {
            [a1, a2, p1, p2] = rk4(a1, a2, p1, p2, dt); // Perform a simulation step
            tail.push(a1, a2); // Record the new angles in history
        },
        /* Create a new pendulum instance with slightly different initial conditions */
        clone: function(conf) {
            if (!conf)
                conf = {};
            let cp2;
            if (p2 === 0.0)
                cp2 = Math.random() * 1e-12;  // Add slight randomness if p2 is zero
            else
                cp2 = p2 * (1 - Math.random() * 1e-10); // Slightly modify p2 if non-zero
            conf.init = [a1, a2, p1, cp2]; // Set initial state for clone
            return new pendulum(conf);
        },
    };
}

function clear(gl) {
    gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
    gl.clear(gl.COLOR_BUFFER_BIT);
}

function draw(gl, webgl, pendulums) {
    let w = gl.canvas.width;
    let h = gl.canvas.height;
    let z = Math.min(w, h);
    let ax = w / z;
    let ay = h / z;
    let d = barLength * 2;
    // Draw tail
    let tail = webgl.tail;
    gl.useProgram(tail.program);
    gl.uniform2f(tail.u_aspect, ax / d, ay / d);
    gl.bindBuffer(gl.ARRAY_BUFFER, webgl.alpha);
    gl.enableVertexAttribArray(tail.a_alpha);
    gl.vertexAttribPointer(tail.a_alpha, 1, gl.FLOAT, false, 0, 0);
    gl.bindBuffer(gl.ARRAY_BUFFER, webgl.tailb);
    gl.enableVertexAttribArray(tail.a_point);
    gl.vertexAttribPointer(tail.a_point, 2, gl.FLOAT, false, 0, 0);
    for (let i = 0; i < pendulums.length; i++) {
        let p = pendulums[i];
        if (p.tail.length) {
            polyline(p.tail, webgl.tailpoly);
            gl.bufferSubData(gl.ARRAY_BUFFER, 0, webgl.tailpoly);
            gl.uniform3fv(tail.u_color, p.tailColor);
            let cutoff = 1 - p.tail.length * 2 / p.tail.v.length;
            gl.uniform1f(tail.u_cutoff, cutoff);
            gl.drawArrays(gl.TRIANGLE_STRIP, 0, p.tail.length * 2);
        }
    }
    // Draw Mass
    let mass = webgl.mass;
    gl.useProgram(mass.program);
    gl.uniform2f(mass.u_aspect, ax, ay);
    gl.bindBuffer(gl.ARRAY_BUFFER, webgl.quad);
    gl.enableVertexAttribArray(mass.a_point);
    gl.vertexAttribPointer(mass.a_point, 2, gl.FLOAT, false, 0, 0);
    for (let i = 0; i < pendulums.length; i++) {
        let p = pendulums[i];
        let [x1, y1, x2, y2] = p.positions();
        x1 *= d / ax;
        y1 *= d / ay;
        x2 *= d / ax;
        y2 *= d / ay;
        gl.uniform3fv(mass.u_color, p.massColor);
        gl.uniform2f(mass.u_center, x1, y1);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
        gl.uniform2f(mass.u_center, x2, y2);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    }
    // Draw Bar
    let bar = webgl.bar;
    gl.useProgram(bar.program);
    gl.uniform2f(bar.u_aspect, ax, ay);
    gl.enableVertexAttribArray(bar.a_point);
    /* Quad buffer still bound from previous draws */
    gl.vertexAttribPointer(bar.a_point, 2, gl.FLOAT, false, 0, 0);
    for (let i = 0; i < pendulums.length; i++) {
        let p = pendulums[i];
        let [x1, y1, x2, y2] = p.positions();
        let [a1, a2, p1, p2] = p.state();
        x1 *= d / ax;
        y1 *= d / ay;
        gl.uniform3fv(bar.u_color, p.massColor);
        gl.uniform2f(bar.u_attach, 0, 0);
        gl.uniform1f(bar.u_angle, a1 - Math.PI / 2);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
        gl.uniform2f(bar.u_attach, x1, y1);
        gl.uniform1f(bar.u_angle, a2 - Math.PI / 2);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    }
};

function glRenderer(gl, tailLen) {
    let webgl = {};
    gl.clearColor(1, 1, 1, 1);
    gl.enable(gl.BLEND); // Enable blending for transparency
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); // Set blend function for transparency

    // Create and initialize a buffer for a quad
    webgl.quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, webgl.quad);
    gl.bufferData(gl.ARRAY_BUFFER, quad, gl.STATIC_DRAW);

    // Create and initialize a buffer for the tail
    webgl.tailb = gl.createBuffer();
    webgl.tailpoly = new Float32Array(tailLen * 4); // Allocate memory for the tail vertices
    gl.bindBuffer(gl.ARRAY_BUFFER, webgl.tailb);
    gl.bufferData(gl.ARRAY_BUFFER, webgl.tailpoly.byteLength, gl.STREAM_DRAW); // Use STREAM_DRAW for dynamic data

    // Create and initialize a buffer for the tail alpha values
    webgl.alpha = gl.createBuffer();
    let alpha = new Float32Array(tailLen * 2);
    for (let i = 0; i < alpha.length; i++) {
        let v = (i + 1) / alpha.length;
        alpha[i] = 1 - v;
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, webgl.alpha);
    gl.bufferData(gl.ARRAY_BUFFER, alpha, gl.STATIC_DRAW);

    // Compile shaders for mass, bar, and tail
    webgl.mass = compile(gl, massShader.vert, massShader.frag);
    webgl.bar  = compile(gl, barShader.vert, barShader.frag);
    webgl.tail = compile(gl, tailShader.vert, tailShader.frag);

    webgl.renderAll = function(pendulums) {
        clear(gl); // Clear the screen
        draw(gl, webgl, pendulums); // Draw all pendulums
    };
    return webgl;
}

(function() {
    let state = [new pendulum()];
    let params = new URL(document.location);
    let useWebGL = params.searchParams.get("webgl") !== '0';
    let canvas = document.getElementById('canvas');
    let running = true;
    let gl = useWebGL ? canvas.getContext('webgl') : null;
    let renderer = null;

    if (!gl) {
        console.log("ERROR: cannot load webgl 3d");
    } else {
        renderer = new glRenderer(gl, tailMax);
    }

    window.addEventListener('keypress', function(e) {
        switch (e.key) {
            case ' ': // SPACE
                running = !running;
                break;
            case 'a': // a
                let color = [Math.random(), Math.random(), Math.random()];
                state.push(new pendulum({tailColor: color}));
                break;
            case 'c': // c
                if (state.length) {
                    let color = [Math.random(), Math.random(), Math.random()];
                    state.push(state[0].clone({tailColor: color}));
                }
                break;
            case 'd': // d
                if (state.length)
                    state.pop();
                break;
        }
    });

    let last = 0.0;
    function cb(t) {
        let dt = Math.min(t - last, dtMax);
        let ww = window.innerWidth;
        let wh = window.innerHeight;
        if (canvas.width != ww || canvas.height != wh) {
            /* Only resize when necessary */
            canvas.width = ww;
            canvas.height = wh;
        }
        if (running)
            for (let i = 0; i < state.length; i++)
                state[i].step(dt / 1000.0);
        clear(gl);
        renderer.renderAll(state);
        last = t;
        window.requestAnimationFrame(cb);
    }

    window.requestAnimationFrame(cb);
}());
