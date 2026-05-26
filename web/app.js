/**
 * Antigravity Graph Drawing & Analysis - Main Application Script
 * Core Canvas Engine, Regression Math Solvers, and UI Controller
 */

// ==========================================================================
// 1. MATHEMATICAL REGRESSION & SOLVERS ENGINE
// ==========================================================================
const MathEngine = {
  /**
   * Calculates linear regression (y = mx + c)
   * @param {Array} points Array of {x, y} points
   * @returns {Object|null} {m, c, r2, formula} or null
   */
  linearRegression(points) {
    const N = points.length;
    if (N < 2) return null;

    let sumX = 0, sumY = 0, sumXX = 0, sumYY = 0, sumXY = 0;
    for (let i = 0; i < N; i++) {
      const p = points[i];
      sumX += p.x;
      sumY += p.y;
      sumXX += p.x * p.x;
      sumYY += p.y * p.y;
      sumXY += p.x * p.y;
    }

    const denominator = N * sumXX - sumX * sumX;
    // Handle vertical line / collinear X values
    if (Math.abs(denominator) < 1e-12) return null;

    const m = (N * sumXY - sumX * sumY) / denominator;
    const c = (sumY - m * sumX) / N;

    // R2 Coefficient
    const yMean = sumY / N;
    let ssTot = 0;
    let ssRes = 0;
    for (let i = 0; i < N; i++) {
      const p = points[i];
      const yFit = m * p.x + c;
      ssTot += Math.pow(p.y - yMean, 2);
      ssRes += Math.pow(p.y - yFit, 2);
    }

    const r2 = ssTot === 0 ? 1.0 : 1 - (ssRes / ssTot);

    // Format Equation
    const mStr = this.formatNumber(m);
    const cSign = c >= 0 ? '+' : '-';
    const cStr = this.formatNumber(Math.abs(c));
    const formula = `y = ${mStr}x ${cSign} ${cStr}`;

    return { m, c, r2, formula };
  },

  /**
   * Calculates quadratic regression (y = ax^2 + bx + c) using Cramer's Rule for 3x3 matrix
   * @param {Array} points Array of {x, y} points
   * @returns {Object|null} {a, b, c, r2, formula} or null
   */
  quadraticRegression(points) {
    const N = points.length;
    if (N < 3) return null;

    let sumX = 0, sumY = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0;
    let sumXY = 0, sumX2Y = 0;

    for (let i = 0; i < N; i++) {
      const p = points[i];
      const x = p.x;
      const y = p.y;
      const x2 = x * x;

      sumX += x;
      sumY += y;
      sumX2 += x2;
      sumX3 += x2 * x;
      sumX4 += x2 * x2;
      sumXY += x * y;
      sumX2Y += x2 * y;
    }

    // Set up matrix system M * A = B
    // [ sumX4  sumX3  sumX2 ] [ a ]   [ sumX2Y ]
    // [ sumX3  sumX2  sumX  ] [ b ] = [ sumXY  ]
    // [ sumX2  sumX   N     ] [ c ]   [ sumY   ]
    const M = [
      [sumX4, sumX3, sumX2],
      [sumX3, sumX2, sumX],
      [sumX2, sumX, N]
    ];
    const B = [sumX2Y, sumXY, sumY];

    // 3x3 Determinant Helper
    const det3x3 = (matrix) => {
      return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
             matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
             matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    };

    const detM = det3x3(M);
    if (Math.abs(detM) < 1e-12) return null; // No unique quadratic fit possible

    // Cramer's rule matrices
    const Ma = [
      [B[0], M[0][1], M[0][2]],
      [B[1], M[1][1], M[1][2]],
      [B[2], M[2][1], M[2][2]]
    ];
    const Mb = [
      [M[0][0], B[0], M[0][2]],
      [M[1][0], B[1], M[1][2]],
      [M[2][0], B[2], M[2][2]]
    ];
    const Mc = [
      [M[0][0], M[0][1], B[0]],
      [M[1][0], M[1][1], B[1]],
      [M[2][0], M[2][1], B[2]]
    ];

    const a = det3x3(Ma) / detM;
    const b = det3x3(Mb) / detM;
    const c = det3x3(Mc) / detM;

    // R2 Coefficient
    const yMean = sumY / N;
    let ssTot = 0;
    let ssRes = 0;
    for (let i = 0; i < N; i++) {
      const p = points[i];
      const yFit = a * p.x * p.x + b * p.x + c;
      ssTot += Math.pow(p.y - yMean, 2);
      ssRes += Math.pow(p.y - yFit, 2);
    }

    const r2 = ssTot === 0 ? 1.0 : 1 - (ssRes / ssTot);

    // Format equation string
    const aStr = this.formatNumber(a);
    const bSign = b >= 0 ? '+' : '-';
    const bStr = this.formatNumber(Math.abs(b));
    const cSign = c >= 0 ? '+' : '-';
    const cStr = this.formatNumber(Math.abs(c));
    const formula = `y = ${aStr}x² ${bSign} ${bStr}x ${cSign} ${cStr}`;

    return { a, b, c, r2, formula };
  },

  /**
   * Computes step-by-step gradients between consecutive sorted points
   * @param {Array} points List of points
   * @returns {Array} List of {p1, p2, dx, dy, m} gradient info
   */
  calculateStepGradients(points) {
    if (points.length < 2) return [];

    // Sort a copy of points by X coordinate
    const sorted = [...points].sort((p1, p2) => p1.x - p2.x);
    const gradients = [];

    for (let i = 0; i < sorted.length - 1; i++) {
      const p1 = sorted[i];
      const p2 = sorted[i + 1];
      const dx = p2.x - p1.x;
      const dy = p2.y - p1.y;
      
      let m = NaN;
      if (Math.abs(dx) > 1e-12) {
        m = dy / dx;
      }

      gradients.push({
        p1,
        p2,
        dx,
        dy,
        m
      });
    }

    return gradients;
  },

  /**
   * Formats floating-point numbers beautifully
   */
  formatNumber(val, decimals = 3) {
    if (val === undefined || isNaN(val)) return 'NaN';
    if (Math.abs(val) < 1e-10) return '0';
    if (Math.abs(val) >= 1e6 || Math.abs(val) < 1e-3) {
      return val.toExponential(decimals);
    }
    return Number(val.toFixed(decimals)).toString();
  }
};


// ==========================================================================
// 2. INTERACTIVE GRAPH CANVAS SYSTEM
// ==========================================================================
class GraphCanvas {
  constructor(canvasElement, onHoverCallback, onDoubleClickCallback) {
    this.canvas = canvasElement;
    this.ctx = canvasElement.getContext('2d');
    this.onHover = onHoverCallback;
    this.onDoubleClick = onDoubleClickCallback;

    // Viewport transformations (defaults)
    this.scale = 40;        // pixels per math unit
    this.offsetX = 0;       // X displacement of center from origin
    this.offsetY = 0;       // Y displacement of center from origin

    // Drag/Pan states
    this.isPanning = false;
    this.startX = 0;
    this.startY = 0;

    // Hover state
    this.hoveredPointIndex = null;
    this.mouseMathX = 0;
    this.mouseMathY = 0;

    // Bind event listeners
    this.initEvents();
    this.resize();
  }

  resize() {
    const parent = this.canvas.parentElement;
    this.canvas.width = parent.clientWidth;
    this.canvas.height = parent.clientHeight;
    this.draw();
  }

  initEvents() {
    // Resize Listener
    window.addEventListener('resize', () => this.resize());

    // Mouse Panning (Drag)
    this.canvas.addEventListener('mousedown', (e) => {
      if (e.button === 0) { // Left click
        this.isPanning = true;
        this.startX = e.clientX - this.offsetX;
        this.startY = e.clientY - this.offsetY;
      }
    });

    window.addEventListener('mousemove', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const screenX = e.clientX - rect.left;
      const screenY = e.clientY - rect.top;

      // Update current mouse math coordinates for HUD
      const math = this.screenToMath(screenX, screenY);
      this.mouseMathX = math.x;
      this.mouseMathY = math.y;

      if (this.isPanning) {
        this.offsetX = e.clientX - this.startX;
        this.offsetY = e.clientY - this.startY;
        this.draw();
      } else {
        // Trigger hover checks
        this.checkHover(screenX, screenY);
      }

      // Fire coordinate update
      if (this.onHover) this.onHover(this.mouseMathX, this.mouseMathY);
    });

    window.addEventListener('mouseup', () => {
      this.isPanning = false;
    });

    // Zooming via Scroll Wheel (Centered at mouse position)
    this.canvas.addEventListener('wheel', (e) => {
      e.preventDefault();
      const rect = this.canvas.getBoundingClientRect();
      const screenX = e.clientX - rect.left;
      const screenY = e.clientY - rect.top;

      const zoomFactor = e.deltaY < 0 ? 1.15 : 0.85;
      this.zoomAt(screenX, screenY, zoomFactor);
    }, { passive: false });

    // Double-click to Edit / Add Point
    this.canvas.addEventListener('dblclick', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const screenX = e.clientX - rect.left;
      const screenY = e.clientY - rect.top;

      const clickMath = this.screenToMath(screenX, screenY);
      
      // Let app controller handle double click at coordinate clickMath
      if (this.onDoubleClick) {
        this.onDoubleClick(clickMath, this.hoveredPointIndex);
      }
    });
  }

  /**
   * Performs a zoom focused on a screen coordinate (e.g. mouse cursor)
   */
  zoomAt(screenX, screenY, factor) {
    const minScale = 1.5;
    const maxScale = 5000;
    const newScale = Math.max(minScale, Math.min(maxScale, this.scale * factor));

    if (newScale === this.scale) return;

    // Preserve math coordinate position under mouse cursor
    const math = this.screenToMath(screenX, screenY);
    
    this.scale = newScale;
    this.offsetX = screenX - this.canvas.width / 2 - math.x * this.scale;
    this.offsetY = screenY - this.canvas.height / 2 + math.y * this.scale;

    this.draw();
  }

  /**
   * Recenter viewport to origin (0, 0)
   */
  recenter() {
    this.scale = 40;
    this.offsetX = 0;
    this.offsetY = 0;
    this.draw();
  }

  /**
   * Translates screen pixels to math coordinates
   */
  screenToMath(sx, sy) {
    const cx = this.canvas.width / 2;
    const cy = this.canvas.height / 2;
    return {
      x: (sx - cx - this.offsetX) / this.scale,
      y: (cy + this.offsetY - sy) / this.scale
    };
  }

  /**
   * Translates math coordinates to screen pixels
   */
  mathToScreen(mx, my) {
    const cx = this.canvas.width / 2;
    const cy = this.canvas.height / 2;
    return {
      x: cx + this.offsetX + mx * this.scale,
      y: cy + this.offsetY - my * this.scale
    };
  }

  /**
   * Scans if mouse is hovering near a plotted coordinate
   */
  checkHover(sx, sy) {
    if (!window.app || !window.app.points.length) {
      this.hoveredPointIndex = null;
      return;
    }

    const hoverThreshold = 12; // pixels
    let foundIndex = null;

    for (let i = 0; i < window.app.points.length; i++) {
      const p = window.app.points[i];
      const screenPos = this.mathToScreen(p.x, p.y);
      const dist = Math.hypot(sx - screenPos.x, sy - screenPos.y);

      if (dist < hoverThreshold) {
        foundIndex = i;
        break;
      }
    }

    if (this.hoveredPointIndex !== foundIndex) {
      this.hoveredPointIndex = foundIndex;
      this.draw();
    }
  }

  /**
   * Renders the dynamic grid, axes, labels, fits, and points
   */
  draw() {
    // 1. Clear background
    this.ctx.fillStyle = '#07090e';
    this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

    // 2. Draw Dynamic Grid Lines
    this.drawGrid();

    // 3. Draw Main Axes (X & Y axes)
    this.drawAxes();

    // 4. Draw Regression Fits
    if (window.app) {
      if (window.app.linearFitActive) this.drawLinearFit();
      if (window.app.quadraticFitActive) this.drawQuadraticFit();
      this.drawPoints();
    }
  }

  /**
   * Calculates ideal grid subdivision step dynamically based on scale
   */
  getGridSpacing() {
    const targetPixelSpacing = 80; // Ideal distance between lines
    const targetMathSpacing = targetPixelSpacing / this.scale;
    
    // Find power of 10 just below target spacing
    const powerOf10 = Math.pow(10, Math.floor(Math.log10(targetMathSpacing)));
    const ratio = targetMathSpacing / powerOf10;

    if (ratio < 2) return powerOf10;
    if (ratio < 5) return 2 * powerOf10;
    return 5 * powerOf10;
  }

  drawGrid() {
    const step = this.getGridSpacing();
    const cx = this.canvas.width / 2;
    const cy = this.canvas.height / 2;

    // Calculate bounding box in math units
    const topLeft = this.screenToMath(0, 0);
    const bottomRight = this.screenToMath(this.canvas.width, this.canvas.height);

    // Round bounds to nearest step intervals
    const startX = Math.floor(topLeft.x / step) * step;
    const endX = Math.ceil(bottomRight.x / step) * step;
    const startY = Math.floor(bottomRight.y / step) * step;
    const endY = Math.ceil(topLeft.y / step) * step;

    // Format decimals dynamically to fit grid density
    const decimals = Math.max(0, -Math.floor(Math.log10(step)));

    this.ctx.lineWidth = 1;
    this.ctx.font = '10px "JetBrains Mono", monospace';
    this.ctx.textBaseline = 'middle';

    // Draw Vertical Grid Lines
    for (let x = startX; x <= endX; x += step) {
      if (Math.abs(x) < 1e-12) continue; // Skip primary Y-axis

      const screenPos = this.mathToScreen(x, 0);
      
      // Draw grid line
      this.ctx.strokeStyle = 'rgba(255, 255, 255, 0.04)';
      this.ctx.beginPath();
      this.ctx.moveTo(screenPos.x, 0);
      this.ctx.lineTo(screenPos.x, this.canvas.height);
      this.ctx.stroke();

      // Axis label position (stick to bottom/top or axis)
      let labelY = cy + this.offsetY + 12;
      this.ctx.textAlign = 'center';
      this.ctx.fillStyle = 'rgba(255, 255, 255, 0.35)';

      // Draw label coordinate text
      if (labelY < 12) labelY = 12;
      if (labelY > this.canvas.height - 18) labelY = this.canvas.height - 18;
      this.ctx.fillText(x.toFixed(decimals), screenPos.x, labelY);
    }

    // Draw Horizontal Grid Lines
    for (let y = startY; y <= endY; y += step) {
      if (Math.abs(y) < 1e-12) continue; // Skip primary X-axis

      const screenPos = this.mathToScreen(0, y);

      // Draw grid line
      this.ctx.strokeStyle = 'rgba(255, 255, 255, 0.04)';
      this.ctx.beginPath();
      this.ctx.moveTo(0, screenPos.y);
      this.ctx.lineTo(this.canvas.width, screenPos.y);
      this.ctx.stroke();

      // Axis label position (stick to side or axis)
      let labelX = cx + this.offsetX - 12;
      this.ctx.textAlign = 'right';
      this.ctx.fillStyle = 'rgba(255, 255, 255, 0.35)';

      if (labelX < 12 + 10) labelX = 12 + 10;
      if (labelX > this.canvas.width - 20) labelX = this.canvas.width - 20;
      this.ctx.fillText(y.toFixed(decimals), labelX, screenPos.y);
    }
  }

  drawAxes() {
    const cx = this.canvas.width / 2;
    const cy = this.canvas.height / 2;
    const axX = cx + this.offsetX;
    const axY = cy + this.offsetY;

    this.ctx.strokeStyle = 'rgba(255, 255, 255, 0.12)';
    this.ctx.lineWidth = 2;

    // X axis
    if (axY >= 0 && axY <= this.canvas.height) {
      this.ctx.beginPath();
      this.ctx.moveTo(0, axY);
      this.ctx.lineTo(this.canvas.width, axY);
      this.ctx.stroke();
    }

    // Y axis
    if (axX >= 0 && axX <= this.canvas.width) {
      this.ctx.beginPath();
      this.ctx.moveTo(axX, 0);
      this.ctx.lineTo(axX, this.canvas.height);
      this.ctx.stroke();
    }

    // Origin indicator marker label
    if (axX >= 0 && axX <= this.canvas.width && axY >= 0 && axY <= this.canvas.height) {
      this.ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
      this.ctx.font = '10px "JetBrains Mono", monospace';
      this.ctx.textAlign = 'right';
      this.ctx.textBaseline = 'top';
      this.ctx.fillText('0', axX - 6, axY + 6);
    }
  }

  drawLinearFit() {
    const fit = MathEngine.linearRegression(window.app.points);
    if (!fit) return;

    const boundsLeft = this.screenToMath(0, 0);
    const boundsRight = this.screenToMath(this.canvas.width, this.canvas.height);

    const startX = boundsLeft.x;
    const endX = boundsRight.x;

    const startY = fit.m * startX + fit.c;
    const endY = fit.m * endX + fit.c;

    const screenStart = this.mathToScreen(startX, startY);
    const screenEnd = this.mathToScreen(endX, endY);

    // Glowing style for fit line
    this.ctx.save();
    this.ctx.strokeStyle = 'rgba(59, 130, 246, 0.3)';
    this.ctx.lineWidth = 5;
    this.ctx.beginPath();
    this.ctx.moveTo(screenStart.x, screenStart.y);
    this.ctx.lineTo(screenEnd.x, screenEnd.y);
    this.ctx.stroke();

    // Solid core line
    this.ctx.strokeStyle = '#3b82f6';
    this.ctx.lineWidth = 2.5;
    this.ctx.stroke();
    this.ctx.restore();
  }

  drawQuadraticFit() {
    const fit = MathEngine.quadraticRegression(window.app.points);
    if (!fit) return;

    const steps = 150; // granularity
    const dx = this.canvas.width / steps;

    this.ctx.save();
    this.ctx.strokeStyle = 'rgba(168, 85, 247, 0.28)';
    this.ctx.lineWidth = 5;
    this.ctx.beginPath();

    // Trace curve screen point by screen point
    for (let i = 0; i <= steps; i++) {
      const sx = i * dx;
      const math = this.screenToMath(sx, 0);
      const my = fit.a * math.x * math.x + fit.b * math.x + fit.c;
      const screenPos = this.mathToScreen(math.x, my);

      if (i === 0) {
        this.ctx.moveTo(screenPos.x, screenPos.y);
      } else {
        this.ctx.lineTo(screenPos.x, screenPos.y);
      }
    }
    this.ctx.stroke();

    this.ctx.strokeStyle = '#a855f7';
    this.ctx.lineWidth = 2.5;
    this.ctx.stroke();
    this.ctx.restore();
  }

  drawPoints() {
    for (let i = 0; i < window.app.points.length; i++) {
      const p = window.app.points[i];
      const screen = this.mathToScreen(p.x, p.y);
      const isHovered = (this.hoveredPointIndex === i);

      this.ctx.save();

      // Glowing outer ring for point
      const glowRadius = isHovered ? 12 : 6;
      const pointRadius = isHovered ? 6.5 : 4.5;
      const alpha = isHovered ? 0.45 : 0.25;

      const gradient = this.ctx.createRadialGradient(
        screen.x, screen.y, 1, 
        screen.x, screen.y, glowRadius
      );
      gradient.addColorStop(0, `rgba(0, 242, 254, ${alpha * 2})`);
      gradient.addColorStop(0.3, `rgba(0, 242, 254, ${alpha})`);
      gradient.addColorStop(1, 'rgba(0, 242, 254, 0)');

      this.ctx.fillStyle = gradient;
      this.ctx.beginPath();
      this.ctx.arc(screen.x, screen.y, glowRadius, 0, 2 * Math.PI);
      this.ctx.fill();

      // Solid central core dot
      this.ctx.fillStyle = '#00f2fe';
      this.ctx.beginPath();
      this.ctx.arc(screen.x, screen.y, pointRadius, 0, 2 * Math.PI);
      this.ctx.fill();

      this.ctx.strokeStyle = '#ffffff';
      this.ctx.lineWidth = 1.5;
      this.ctx.stroke();

      // Dynamic tooltip overlay if hovered
      if (isHovered) {
        const text = `(${MathEngine.formatNumber(p.x)}, ${MathEngine.formatNumber(p.y)})`;
        this.ctx.font = '700 11px "JetBrains Mono", monospace';
        this.ctx.fillStyle = '#ffffff';
        this.ctx.textAlign = 'left';
        this.ctx.textBaseline = 'bottom';
        
        // Draw coordinate text hovering above the dot
        this.ctx.shadowColor = 'rgba(0, 0, 0, 0.8)';
        this.ctx.shadowBlur = 4;
        this.ctx.fillText(text, screen.x + 10, screen.y - 10);
      }

      this.ctx.restore();
    }
  }
}


// ==========================================================================
// 3. APPLICATION CONTROLLER (STATE & EVENT BINDINGS)
// ==========================================================================
class AppController {
  constructor() {
    this.points = [];
    this.linearFitActive = false;
    this.quadraticFitActive = false;

    // Fixed Interval Settings
    this.fixedIntervalActive = false;
    this.fixedXStart = 1;
    this.fixedXInterval = 1;

    this.initDOMElements();
    this.initCanvas();
    this.bindEvents();
    this.updateUI();
  }

  initDOMElements() {
    // Buttons
    this.btnOpenAdd = document.getElementById('btn-open-add-dialog');
    this.btnCloseModal = document.getElementById('btn-close-modal');
    this.btnCancelModal = document.getElementById('btn-cancel-modal');
    this.btnSubmitPoint = document.getElementById('btn-submit-point');
    this.btnClearAll = document.getElementById('btn-clear-all');
    this.btnResetView = document.getElementById('btn-reset-view');

    this.btnToggleLinear = document.getElementById('btn-toggle-linear');
    this.btnToggleQuadratic = document.getElementById('btn-toggle-quadratic');

    // Float Canvas Buttons
    this.btnZoomIn = document.getElementById('btn-zoom-in');
    this.btnZoomOut = document.getElementById('btn-zoom-out');
    this.btnRecenter = document.getElementById('btn-recenter');

    // HUD overlays
    this.coordHUD = document.getElementById('coord-hud');

    // Interval DOMs
    this.intervalToggle = document.getElementById('interval-toggle');
    this.intervalPanel = document.getElementById('interval-inputs-panel');
    this.inputIntervalStart = document.getElementById('interval-start');
    this.inputIntervalStep = document.getElementById('interval-step');

    // Modal DOMs
    this.modal = document.getElementById('modal-container');
    this.modalTitle = document.getElementById('modal-title');
    this.formPoint = document.getElementById('point-form');
    this.inputEditIndex = document.getElementById('edit-index');
    this.inputX = document.getElementById('modal-input-x');
    this.inputY = document.getElementById('modal-input-y');
    this.modalXGroup = document.getElementById('modal-x-group');
    this.errorX = document.getElementById('error-x');
    this.errorY = document.getElementById('error-y');

    // Sidebar lists & Readouts
    this.pointsList = document.getElementById('points-list');
    this.pointsEmptyState = document.getElementById('points-empty-state');
    this.pointsCounter = document.getElementById('points-counter');

    this.fitReadoutPanel = document.getElementById('fit-readout-panel');
    this.linearReadout = document.getElementById('linear-readout');
    this.linearFormula = document.getElementById('linear-formula');
    this.linearR2 = document.getElementById('linear-r2');

    this.quadraticReadout = document.getElementById('quadratic-readout');
    this.quadraticFormula = document.getElementById('quadratic-formula');
    this.quadraticR2 = document.getElementById('quadratic-r2');

    this.gradientOverall = document.getElementById('gradient-overall');
    this.stepGradientsList = document.getElementById('step-gradients-list');
  }

  initCanvas() {
    const canvasElement = document.getElementById('graph-canvas');
    this.canvasController = new GraphCanvas(
      canvasElement,
      // OnHover coordinate updater
      (mx, my) => {
        this.coordHUD.textContent = `X: ${MathEngine.formatNumber(mx, 4)} | Y: ${MathEngine.formatNumber(my, 4)}`;
      },
      // OnDoubleClick handler
      (clickMath, hoveredPointIndex) => {
        if (hoveredPointIndex !== null) {
          // Edit hovered point
          this.openCoordinateModal(hoveredPointIndex);
        } else {
          // Add new point at double click location
          this.openCoordinateModal(-1, clickMath.x, clickMath.y);
        }
      }
    );
  }

  bindEvents() {
    // Add point triggers
    this.btnOpenAdd.addEventListener('click', () => this.openCoordinateModal(-1));
    this.btnCloseModal.addEventListener('click', () => this.closeCoordinateModal());
    this.btnCancelModal.addEventListener('click', () => this.closeCoordinateModal());
    this.formPoint.addEventListener('submit', (e) => {
      e.preventDefault();
      this.handleSubmitPoint();
    });

    // Sidebar actions
    this.btnClearAll.addEventListener('click', () => this.clearAllPoints());
    this.btnResetView.addEventListener('click', () => this.canvasController.recenter());

    // Float actions
    this.btnZoomIn.addEventListener('click', () => {
      const cx = this.canvasController.canvas.width / 2;
      const cy = this.canvasController.canvas.height / 2;
      this.canvasController.zoomAt(cx, cy, 1.25);
    });
    this.btnZoomOut.addEventListener('click', () => {
      const cx = this.canvasController.canvas.width / 2;
      const cy = this.canvasController.canvas.height / 2;
      this.canvasController.zoomAt(cx, cy, 0.8);
    });
    this.btnRecenter.addEventListener('click', () => this.canvasController.recenter());

    // Switch/Modes Actions
    this.intervalToggle.addEventListener('change', (e) => {
      this.fixedIntervalActive = e.target.checked;
      if (this.fixedIntervalActive) {
        this.intervalPanel.classList.remove('hidden');
        this.recalculateFixedIntervals();
      } else {
        this.intervalPanel.classList.add('hidden');
      }
      this.updateUI();
    });

    this.inputIntervalStart.addEventListener('input', () => {
      this.fixedXStart = parseFloat(this.inputIntervalStart.value) || 1;
      if (this.fixedIntervalActive) {
        this.recalculateFixedIntervals();
        this.updateUI();
      }
    });

    this.inputIntervalStep.addEventListener('input', () => {
      this.fixedXInterval = Math.max(0.0001, parseFloat(this.inputIntervalStep.value) || 1);
      if (this.fixedIntervalActive) {
        this.recalculateFixedIntervals();
        this.updateUI();
      }
    });

    // Toggle regression fits
    this.btnToggleLinear.addEventListener('click', () => {
      this.linearFitActive = !this.linearFitActive;
      this.btnToggleLinear.setAttribute('data-active', this.linearFitActive);
      this.updateUI();
    });

    this.btnToggleQuadratic.addEventListener('click', () => {
      this.quadraticFitActive = !this.quadraticFitActive;
      this.btnToggleQuadratic.setAttribute('data-active', this.quadraticFitActive);
      this.updateUI();
    });
  }

  /**
   * Opens the coordinates input glass modal
   * @param {number} editIndex index of point to edit, or -1 to add a new point
   * @param {number} [prefX] preset X coordinate (from click)
   * @param {number} [prefY] preset Y coordinate (from click)
   */
  openCoordinateModal(editIndex, prefX, prefY) {
    this.errorX.textContent = '';
    this.errorY.textContent = '';
    this.inputEditIndex.value = editIndex;

    if (editIndex >= 0) {
      // EDIT MODE
      const point = this.points[editIndex];
      this.modalTitle.textContent = `Edit Point #${editIndex + 1}`;
      this.inputX.value = point.x;
      this.inputY.value = point.y;
      this.btnSubmitPoint.textContent = 'Save Changes';

      // Hide or lock X group in Fixed Interval Mode
      if (this.fixedIntervalActive) {
        this.modalXGroup.style.display = 'none';
        this.inputX.removeAttribute('required');
      } else {
        this.modalXGroup.style.display = 'flex';
        this.inputX.setAttribute('required', 'true');
      }
    } else {
      // ADD MODE
      this.modalTitle.textContent = 'Add Point';
      this.btnSubmitPoint.textContent = 'Plot Point';

      if (this.fixedIntervalActive) {
        // Locked X spacing
        const nextX = this.fixedXStart + this.points.length * this.fixedXInterval;
        this.modalXGroup.style.display = 'none';
        this.inputX.value = nextX;
        this.inputX.removeAttribute('required');
      } else {
        this.modalXGroup.style.display = 'flex';
        this.inputX.setAttribute('required', 'true');
        this.inputX.value = prefX !== undefined ? MathEngine.formatNumber(prefX, 3) : '';
      }

      this.inputY.value = prefY !== undefined ? MathEngine.formatNumber(prefY, 3) : '';
    }

    this.modal.classList.add('active');
    // Auto focus on Y input since it's always required and X is often automatic or filled
    setTimeout(() => this.inputY.focus(), 150);
  }

  closeCoordinateModal() {
    this.modal.classList.remove('active');
  }

  /**
   * Handles coordinate form submit (Adds new point or updates edited one)
   */
  handleSubmitPoint() {
    const editIndex = parseInt(this.inputEditIndex.value);
    const yVal = parseFloat(this.inputY.value);
    
    if (isNaN(yVal)) {
      this.errorY.textContent = 'Please enter a valid numeric value.';
      return;
    }

    let xVal;
    if (this.fixedIntervalActive) {
      if (editIndex >= 0) {
        xVal = this.points[editIndex].x; // X doesn't change for edit in interval mode
      } else {
        xVal = this.fixedXStart + this.points.length * this.fixedXInterval;
      }
    } else {
      xVal = parseFloat(this.inputX.value);
      if (isNaN(xVal)) {
        this.errorX.textContent = 'Please enter a valid numeric X coordinate.';
        return;
      }
    }

    if (editIndex >= 0) {
      // Update existing
      this.points[editIndex] = { x: xVal, y: yVal };
    } else {
      // Add new
      this.points.push({ x: xVal, y: yVal });
    }

    this.closeCoordinateModal();
    this.updateUI();
  }

  /**
   * Recalculates X positions of points in Fixed Interval Mode
   */
  recalculateFixedIntervals() {
    for (let i = 0; i < this.points.length; i++) {
      this.points[i].x = this.fixedXStart + i * this.fixedXInterval;
    }
  }

  deletePoint(index) {
    this.points.splice(index, 1);
    if (this.fixedIntervalActive) {
      // Re-space subsequent X coordinates
      this.recalculateFixedIntervals();
    }
    this.updateUI();
  }

  clearAllPoints() {
    this.points = [];
    this.updateUI();
  }

  /**
   * Refreshes the side list, statistics, equations readouts, and redraws the canvas
   */
  updateUI() {
    const N = this.points.length;

    // 1. Update points list header
    this.pointsCounter.textContent = `${N} Point${N !== 1 ? 's' : ''}`;

    // 2. Render Plotted Points List
    this.pointsList.innerHTML = '';
    if (N === 0) {
      this.pointsList.appendChild(this.pointsEmptyState);
      this.pointsEmptyState.style.display = 'flex';
    } else {
      this.pointsEmptyState.style.display = 'none';

      for (let i = 0; i < N; i++) {
        const p = this.points[i];
        
        const item = document.createElement('div');
        item.className = 'point-item';
        item.role = 'listitem';

        const coordSpan = document.createElement('span');
        coordSpan.className = 'point-coord';
        
        if (this.fixedIntervalActive) {
          coordSpan.innerHTML = `<span class="coord-label">X[${i+1}]</span>${MathEngine.formatNumber(p.x, 2)} → <span class="coord-label">Y</span>${MathEngine.formatNumber(p.y, 4)}`;
        } else {
          coordSpan.innerHTML = `<span class="coord-label">P${i+1}:</span>(${MathEngine.formatNumber(p.x, 3)}, ${MathEngine.formatNumber(p.y, 3)})`;
        }

        const actions = document.createElement('div');
        actions.className = 'point-item-actions';

        // Edit button
        const btnEdit = document.createElement('button');
        btnEdit.className = 'btn-point-action btn-edit';
        btnEdit.title = 'Edit values';
        btnEdit.innerHTML = `
          <svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" stroke-width="2.5">
            <path d="M11 4H4a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h14a2 2 0 0 0 2-2v-7" stroke-linecap="round"/>
            <path d="M18.5 2.5a2.121 2.121 0 1 1 3 3L12 15l-4 1 1-4 9.5-9.5z" stroke-linecap="round"/>
          </svg>
        `;
        btnEdit.addEventListener('click', (e) => {
          e.stopPropagation();
          this.openCoordinateModal(i);
        });

        // Delete button
        const btnDelete = document.createElement('button');
        btnDelete.className = 'btn-point-action btn-delete';
        btnDelete.title = 'Remove point';
        btnDelete.innerHTML = `
          <svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" stroke-width="2.5">
            <polyline points="3 6 5 6 21 6"/>
            <path d="M19 6v14a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2V6m3 0V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v2"/>
          </svg>
        `;
        btnDelete.addEventListener('click', (e) => {
          e.stopPropagation();
          this.deletePoint(i);
        });

        actions.appendChild(btnEdit);
        actions.appendChild(btnDelete);
        item.appendChild(coordSpan);
        item.appendChild(actions);

        // Edit on clicking list item row itself!
        item.addEventListener('dblclick', () => {
          this.openCoordinateModal(i);
        });

        this.pointsList.appendChild(item);
      }
    }

    // 3. Compute regression and update readout panel
    let showReadouts = false;

    // Linear calculations
    if (this.linearFitActive && N >= 2) {
      const fit = MathEngine.linearRegression(this.points);
      if (fit) {
        showReadouts = true;
        this.linearReadout.style.display = 'block';
        this.linearFormula.textContent = fit.formula;
        this.linearR2.textContent = `R² = ${MathEngine.formatNumber(fit.r2, 4)}`;
        this.gradientOverall.textContent = MathEngine.formatNumber(fit.m, 4);
      } else {
        this.linearReadout.style.display = 'none';
        this.gradientOverall.textContent = 'Undefined (Collinear)';
      }
    } else {
      this.linearReadout.style.display = 'none';
      if (!this.linearFitActive) {
        this.gradientOverall.textContent = 'N/A';
      }
    }

    // Quadratic calculations
    if (this.quadraticFitActive && N >= 3) {
      const fit = MathEngine.quadraticRegression(this.points);
      if (fit) {
        showReadouts = true;
        this.quadraticReadout.style.display = 'block';
        this.quadraticFormula.textContent = fit.formula;
        this.quadraticR2.textContent = `R² = ${MathEngine.formatNumber(fit.r2, 4)}`;
      } else {
        this.quadraticReadout.style.display = 'none';
      }
    } else {
      this.quadraticReadout.style.display = 'none';
    }

    if (showReadouts) {
      this.fitReadoutPanel.classList.remove('hidden');
    } else {
      this.fitReadoutPanel.classList.add('hidden');
    }

    // If linear fit isn't active, show overall gradient calculation by regression if possible
    if (!this.linearFitActive && N >= 2) {
      const fit = MathEngine.linearRegression(this.points);
      this.gradientOverall.textContent = fit ? MathEngine.formatNumber(fit.m, 4) : 'Undefined';
    }

    // 4. Update Consecutive Point Slopes
    this.stepGradientsList.innerHTML = '';
    if (N >= 2) {
      const gradients = MathEngine.calculateStepGradients(this.points);
      for (let i = 0; i < gradients.length; i++) {
        const g = gradients[i];

        const item = document.createElement('div');
        item.className = 'slope-item';

        const pointsSpan = document.createElement('span');
        pointsSpan.className = 'slope-points';
        pointsSpan.textContent = `P${this.points.indexOf(g.p1)+1} ➔ P${this.points.indexOf(g.p2)+1}`;

        const slopeSpan = document.createElement('span');
        if (isNaN(g.m)) {
          slopeSpan.className = 'slope-value slope-zero';
          slopeSpan.textContent = '∞ (Vertical)';
        } else {
          const signClass = g.m > 0 ? 'slope-positive' : (g.m < 0 ? 'slope-negative' : 'slope-zero');
          const prefix = g.m > 0 ? '+' : '';
          slopeSpan.className = `slope-value ${signClass}`;
          slopeSpan.textContent = `${prefix}${MathEngine.formatNumber(g.m, 3)}`;
        }

        item.appendChild(pointsSpan);
        item.appendChild(slopeSpan);
        this.stepGradientsList.appendChild(item);
      }
    } else {
      const emptyText = document.createElement('p');
      emptyText.className = 'empty-text';
      emptyText.textContent = 'Plot at least 2 points to view slopes.';
      this.stepGradientsList.appendChild(emptyText);
    }

    // 5. Redraw the canvas content
    this.canvasController.draw();
  }
}

// Instantiate global app controller once script is fully evaluated
window.addEventListener('DOMContentLoaded', () => {
  window.app = new AppController();
});
