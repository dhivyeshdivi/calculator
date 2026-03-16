'use strict';

/* ═══════════════════════════════════════════════════
   UTILITIES
═══════════════════════════════════════════════════ */
const U = {
  parseFn(expr, vars=['x']) {
    if (!expr||!expr.trim()) throw new Error('Expression cannot be empty.');
    let e = expr.trim()
      .replace(/\^/g,'**')
      .replace(/\bsin\b/g,'Math.sin').replace(/\bcos\b/g,'Math.cos')
      .replace(/\btan\b/g,'Math.tan').replace(/\bexp\b/g,'Math.exp')
      .replace(/\bln\b/g,'Math.log').replace(/\blog\b/g,'Math.log10')
      .replace(/\bsqrt\b/g,'Math.sqrt').replace(/\babs\b/g,'Math.abs')
      .replace(/\bfloor\b/g,'Math.floor').replace(/\bceil\b/g,'Math.ceil')
      .replace(/\bpi\b/gi,'Math.PI')
      .replace(/(?<![a-zA-Z])e(?![a-zA-Z0-9_])/g,'Math.E');
    try { return new Function(...vars, `"use strict";return(${e});`); }
    catch { throw new Error(`Invalid expression: "${expr}". Check syntax.`); }
  },
  safe(fn,...args) {
    const r = fn(...args);
    if (!isFinite(r)) throw new Error(`Non-finite result at (${args.map(v=>v.toPrecision(4)).join(', ')}). Check expression/limits.`);
    return r;
  },
  fmt(n, sig=8) {
    if (!isFinite(n)) return String(n);
    return parseFloat(n.toPrecision(sig)).toString();
  },
  float(id, label) {
    const el=document.getElementById(id), v=parseFloat(el.value);
    if (isNaN(v)) { el.classList.add('has-error'); throw new Error(`"${label}" must be a valid number.`); }
    el.classList.remove('has-error'); return v;
  },
  int(id, label, mn=1, mx=1000) {
    const el=document.getElementById(id), v=parseInt(el.value,10);
    if (isNaN(v)||v<mn||v>mx) { el.classList.add('has-error'); throw new Error(`"${label}" must be integer ${mn}–${mx}.`); }
    el.classList.remove('has-error'); return v;
  },
  str(id, label) {
    const el=document.getElementById(id), v=el.value.trim();
    if (!v) { el.classList.add('has-error'); throw new Error(`"${label}" cannot be empty.`); }
    el.classList.remove('has-error'); return v;
  },
  err(id, msg) {
    const el=document.getElementById(id); if (!el) return;
    const sp=el.querySelector('span')||el;
    sp.textContent=msg; el.classList.add('show');
  },
  clearErr(id) {
    const el=document.getElementById(id); if (!el) return;
    (el.querySelector('span')||el).textContent='';
    el.classList.remove('show');
  },
  show(id) { const el=document.getElementById(id); if(el) el.classList.add('show'); },
  hide(id) { const el=document.getElementById(id); if(el) el.classList.remove('show'); },
  clearInputErr() { document.querySelectorAll('.form-input.has-error').forEach(e=>e.classList.remove('has-error')); }
};

/* ═══════════════════════════════════════════════════
   GAUSS–SEIDEL
═══════════════════════════════════════════════════ */
const GS = {
  f5(n){ return n.toFixed(5); },

  run() {
    U.clearErr('gs-error'); U.clearInputErr();
    try {
      const A = [
        [U.float('a11','a₁₁'), U.float('a12','a₁₂'), U.float('a13','a₁₃')],
        [U.float('a21','a₂₁'), U.float('a22','a₂₂'), U.float('a23','a₂₃')],
        [U.float('a31','a₃₁'), U.float('a32','a₃₂'), U.float('a33','a₃₃')],
      ];
      const b = [U.float('b1','b₁'), U.float('b2','b₂'), U.float('b3','b₃')];
      const x = [U.float('x0','x⁽⁰⁾'), U.float('y0','y⁽⁰⁾'), U.float('z0','z⁽⁰⁾')];

      for(let i=0;i<3;i++) if(A[i][i]===0) throw new Error(`Diagonal element a${i+1}${i+1} is zero — rearrange equations so largest coefficients are on the diagonal.`);

      const TOL=1e-5, MAX_ITER=100;

      const buildFormula = (diagIdx, coeffs, bVal, varNames) => {
        const diag = coeffs[diagIdx];
        let terms = [];
        coeffs.forEach((c, i) => {
          if (i === diagIdx || c === 0) return;
          const abs = Math.abs(c);
          const effectiveSign = c > 0 ? '−' : '+';
          const label = abs === 1 ? varNames[i] : `${abs}${varNames[i]}`;
          terms.push(`${effectiveSign} ${label}`);
        });
        const rhs = terms.length ? `${bVal} ${terms.join(' ')}` : `${bVal}`;
        return `${varNames[diagIdx]} = (1 / ${Math.abs(diag)}) [${rhs}]`;
      };

      const vn = ['x','y','z'];
      const formulaLines = [
        buildFormula(0, A[0], b[0], vn) + '  ... (1)',
        buildFormula(1, A[1], b[1], vn) + '  ... (2)',
        buildFormula(2, A[2], b[2], vn) + '  ... (3)',
      ];
      document.getElementById('gs-formulas').innerHTML =
        formulaLines.map(l => `<div>${l}</div>`).join('');

      const allIters = [];
      for(let iter=1; iter<=MAX_ITER; iter++){
        const prev = [...x];
        const steps = [];

        x[0] = (b[0] - A[0][1]*x[1] - A[0][2]*x[2]) / A[0][0];
        steps.push(this.buildSubLine('x', A[0], b[0], [null, prev[1], prev[2]], 0, x[0]));

        x[1] = (b[1] - A[1][0]*x[0] - A[1][2]*x[2]) / A[1][1];
        steps.push(this.buildSubLine('y', A[1], b[1], [x[0], null, prev[2]], 1, x[1]));

        x[2] = (b[2] - A[2][0]*x[0] - A[2][1]*x[1]) / A[2][2];
        steps.push(this.buildSubLine('z', A[2], b[2], [x[0], x[1], null], 2, x[2]));

        allIters.push({ iter, x: x[0], y: x[1], z: x[2], steps });

        const e1=Math.abs(x[0]-prev[0]), e2=Math.abs(x[1]-prev[1]), e3=Math.abs(x[2]-prev[2]);
        if(e1<TOL && e2<TOL && e3<TOL) break;
      }

      document.getElementById('gs-steps').innerHTML = allIters.map(iter => `
        <div class="iter-block">
          <div class="iter-block-title">Iteration ${iter.iter}:</div>
          ${iter.steps.map(s => `
            <div class="iter-step">
              <div class="iter-step-label">${s.varLabel}</div>
              ${s.lines.map((l,i) => `<div class="iter-step-line ${i===s.lines.length-1?'result':i>0?'sub':''}">${l}</div>`).join('')}
            </div>
          `).join('')}
        </div>
      `).join('');

      document.getElementById('gs-tbody').innerHTML = allIters.map((r,i) => {
        const last = i === allIters.length-1;
        return `<tr>
          <td class="dim">${r.iter}</td>
          <td class="${last?'ca':'hi'}">${this.f5(r.x)}</td>
          <td class="${last?'ca':'hi'}">${this.f5(r.y)}</td>
          <td class="${last?'ca':'hi'}">${this.f5(r.z)}</td>
        </tr>`;
      }).join('');

      document.getElementById('gf1').textContent = this.f5(x[0]);
      document.getElementById('gf2').textContent = this.f5(x[1]);
      document.getElementById('gf3').textContent = this.f5(x[2]);
      document.getElementById('gs-conv-badge').textContent =
        allIters.length < MAX_ITER ? `Converged in ${allIters.length} iteration${allIters.length>1?'s':''}` : `Max ${MAX_ITER} iterations`;

      U.show('gs-results');
      document.getElementById('gs-results').scrollIntoView({behavior:'smooth', block:'start'});
    } catch(e){ U.err('gs-error', e.message); }
  },

  buildSubLine(varName, coeffRow, bVal, vals, skipIdx, result) {
    const diag = Math.abs(coeffRow[skipIdx]);
    const varLabel = `${varName} = (1 / ${diag}) [...]`;
    let terms = [`${bVal}`];
    coeffRow.forEach((c, i) => {
      if (i === skipIdx || c === 0) return;
      const val = vals[i];
      const abs = Math.abs(c);
      const sign = c > 0 ? '−' : '+';
      const label = abs === 1 ? `(${this.f5(val)})` : `${abs}(${this.f5(val)})`;
      terms.push(`${sign} ${label}`);
    });
    const bracket = terms.join(' ');
    const innerVal = coeffRow.reduce((sum, c, i) => {
      if (i === skipIdx) return sum;
      return sum - c * (vals[i] || 0);
    }, bVal);
    return {
      varLabel,
      lines: [
        `${varName} = (1 / ${diag}) [${bracket}]`,
        `   = (1 / ${diag}) [${this.f5(innerVal)}]`,
        `   = ${this.f5(result)}`,
      ]
    };
  },

  reset() {
    ['a11','a12','a13','a21','a22','a23','a31','a32','a33','b1','b2','b3']
      .forEach(id=>{const el=document.getElementById(id);if(el){el.value='';el.classList.remove('has-error');}});
    ['x0','y0','z0'].forEach(id=>{const el=document.getElementById(id);if(el){el.value='0';el.classList.remove('has-error');}});
    U.hide('gs-results'); U.clearErr('gs-error');
    document.getElementById('gs-steps').innerHTML='';
    document.getElementById('gs-formulas').innerHTML='';
  },

  example() {
    Object.entries({
      a11:27, a12:6,  a13:-1, b1:85,
      a21:6,  a22:15, a23:2,  b2:72,
      a31:1,  a32:1,  a33:54, b3:110,
      x0:0, y0:0, z0:0
    }).forEach(([k,v])=>{const el=document.getElementById(k);if(el)el.value=v;});
  }
};

/* ═══════════════════════════════════════════════════
   EULER METHOD
═══════════════════════════════════════════════════ */
const EM = {
  run() {
    U.clearErr('euler-error'); U.clearInputErr();
    try {
      const fExpr=U.str('euler-f','f(x,y)');
      const x0=U.float('euler-x0','x₀'), y0=U.float('euler-y0','y₀');
      const h=U.float('euler-h','h'), n=U.int('euler-n','n',1,500);
      if(h===0) throw new Error('Step size h cannot be zero.');
      const f=U.parseFn(fExpr,['x','y']);
      U.safe(f,x0,y0);

      const rows=[]; let x=x0,y=y0;
      for(let i=0;i<n;i++){
        const fy=U.safe(f,x,y), yn1=y+h*fy;
        rows.push({i,x,y,fy,yn1});
        x=parseFloat((x+h).toPrecision(14)); y=yn1;
      }

      document.getElementById('euler-tbody').innerHTML=rows.map((r,i)=>
        `<tr><td class="dim">${r.i}</td><td class="hi">${U.fmt(r.x)}</td><td>${U.fmt(r.y)}</td><td>${U.fmt(r.fy)}</td><td class="${i===rows.length-1?'ct':''}">${U.fmt(r.yn1)}</td></tr>`
      ).join('');
      const last=rows[rows.length-1];
      document.getElementById('ef-x').textContent=U.fmt(last.x);
      document.getElementById('ef-y').textContent=U.fmt(last.yn1);
      document.getElementById('ef-n').textContent=`${rows.length} step${rows.length!==1?'s':''}`;
      U.show('euler-results');
      document.getElementById('euler-results').scrollIntoView({behavior:'smooth',block:'start'});
    } catch(e){ U.err('euler-error',e.message); }
  },
  reset() {
    ['euler-f','euler-x0','euler-y0','euler-h'].forEach(id=>{const el=document.getElementById(id);if(el){el.value='';el.classList.remove('has-error');}});
    document.getElementById('euler-n').value='5';
    U.hide('euler-results'); U.clearErr('euler-error');
  },
  example() {
    document.getElementById('euler-f').value='x + y';
    document.getElementById('euler-x0').value='0';
    document.getElementById('euler-y0').value='1';
    document.getElementById('euler-h').value='0.1';
    document.getElementById('euler-n').value='5';
  }
};

/* ═══════════════════════════════════════════════════
   TRAPEZOIDAL RULE
═══════════════════════════════════════════════════ */
const TR = {
  run() {
    U.clearErr('trap-error'); U.clearInputErr();
    try {
      const fExpr=U.str('trap-f','f(x)');
      const a=U.float('trap-a','a'), b=U.float('trap-b','b'), n=U.int('trap-n','n',1,2000);
      if(a>=b) throw new Error('Lower limit a must be strictly less than b.');
      const f=U.parseFn(fExpr,['x']);
      const h=(b-a)/n; let wSum=0;
      const rows=[];
      for(let i=0;i<=n;i++){
        const xi=a+i*h, fxi=U.safe(f,xi);
        const coeff=(i===0||i===n)?1:2;
        wSum+=coeff*fxi;
        rows.push({i,xi,fxi,coeff,term:coeff*fxi});
      }
      const I=(h/2)*wSum;
      document.getElementById('trap-h-disp').textContent=U.fmt(h);
      document.getElementById('trap-tbody').innerHTML=rows.map(r=>
        `<tr><td class="dim">${r.i}</td><td class="hi">${U.fmt(r.xi)}</td><td>${U.fmt(r.fxi)}</td><td style="color:${r.coeff===2?'var(--amber)':'var(--rust)'};font-weight:600">${r.coeff}</td><td>${U.fmt(r.term)}</td></tr>`
      ).join('');
      document.getElementById('trap-final').textContent=U.fmt(I);
      U.show('trap-results');
      document.getElementById('trap-results').scrollIntoView({behavior:'smooth',block:'start'});
    } catch(e){ U.err('trap-error',e.message); }
  },
  reset() {
    ['trap-f','trap-a','trap-b'].forEach(id=>{const el=document.getElementById(id);if(el){el.value='';el.classList.remove('has-error');}});
    document.getElementById('trap-n').value='4';
    U.hide('trap-results'); U.clearErr('trap-error');
  },
  example() {
    document.getElementById('trap-f').value='x^2';
    document.getElementById('trap-a').value='0';
    document.getElementById('trap-b').value='1';
    document.getElementById('trap-n').value='4';
  }
};

/* ═══════════════════════════════════════════════════
   DOUBLE INTEGRAL TRAPEZOIDAL
═══════════════════════════════════════════════════ */
const DT = {
  _gridData: null,

  run() {
    U.clearErr('dtrap-error'); U.clearInputErr();
    try {
      const fExpr=U.str('dtrap-f','f(x,y)');
      const a=U.float('dtrap-a','a (lower x)'), b=U.float('dtrap-b','b (upper x)');
      const c=U.float('dtrap-c','c (lower y)'), d=U.float('dtrap-d','d (upper y)');
      const m=U.int('dtrap-m','m (x intervals)',1,100), n=U.int('dtrap-n','n (y intervals)',1,100);
      if(a>=b) throw new Error('x lower limit a must be strictly less than b.');
      if(c>=d) throw new Error('y lower limit c must be strictly less than d.');
      const f=U.parseFn(fExpr,['x','y']);
      U.safe(f,a,c);
      const hx=(b-a)/m, hy=(d-c)/n;
      const xs=[], ys=[];
      for(let i=0;i<=m;i++) xs.push(a+i*hx);
      for(let j=0;j<=n;j++) ys.push(c+j*hy);
      const wx = i => (i===0||i===m) ? 0.5 : 1.0;
      const wy = j => (j===0||j===n) ? 0.5 : 1.0;
      const innerRows=[]; let outerSum=0; const grid=[];
      for(let j=0;j<=n;j++){
        const yj=ys[j], wyj=wy(j); let innerSum=0; const rowF=[];
        for(let i=0;i<=m;i++){ const xi=xs[i], fval=U.safe(f,xi,yj); innerSum+=wx(i)*fval; rowF.push(fval); }
        const Iyj=hx*innerSum; outerSum+=wyj*Iyj;
        innerRows.push({j,yj,wyj,Iyj,contrib:wyj*Iyj}); grid.push(rowF);
      }
      const result=hy*outerSum;
      this._gridData={grid,xs,ys,m,n};
      document.getElementById('dt-hx').textContent=U.fmt(hx);
      document.getElementById('dt-hy').textContent=U.fmt(hy);
      document.getElementById('dt-tbody-y').innerHTML=innerRows.map((r,i)=>{
        const last=i===innerRows.length-1;
        return `<tr><td class="dim">${r.j}</td><td class="hi">${U.fmt(r.yj)}</td><td style="color:var(--violet);font-weight:600">${U.fmt(r.wyj)}</td><td class="${last?'cv':''}">${U.fmt(r.Iyj)}</td><td>${U.fmt(r.contrib)}</td></tr>`;
      }).join('');
      document.getElementById('dt-final').textContent=U.fmt(result);
      document.getElementById('dt-grid-size').textContent=`(${m+1}) × (${n+1}) = ${(m+1)*(n+1)} pts`;
      document.getElementById('dt-toggle-btn').style.display='';
      document.getElementById('dt-grid-wrap').style.display='none';
      document.getElementById('dt-toggle-btn').textContent='Show Full Grid';
      U.show('dtrap-results');
      document.getElementById('dtrap-results').scrollIntoView({behavior:'smooth',block:'start'});
    } catch(e){ U.err('dtrap-error',e.message); }
  },

  toggleGrid() {
    const wrap=document.getElementById('dt-grid-wrap'), btn=document.getElementById('dt-toggle-btn');
    if (wrap.style.display==='none') { this.renderGrid(); wrap.style.display=''; btn.textContent='Hide Full Grid'; }
    else { wrap.style.display='none'; btn.textContent='Show Full Grid'; }
  },

  renderGrid() {
    if (!this._gridData) return;
    const {grid,xs,ys,m,n}=this._gridData;
    const tbl=document.getElementById('dt-grid-table');
    let html='<thead><tr><th>yⱼ \\ xᵢ</th>';
    xs.forEach((x,i)=>{ html+=`<th>x${i}=${U.fmt(x,4)}</th>`; });
    html+='</tr></thead><tbody>';
    grid.forEach((row,j)=>{
      html+=`<tr><td class="hi" style="background:var(--ink);color:var(--cream2)">y${j}=${U.fmt(ys[j],4)}</td>`;
      row.forEach(v=>{ html+=`<td>${U.fmt(v,5)}</td>`; });
      html+='</tr>';
    });
    tbl.innerHTML=html+'</tbody>';
  },

  reset() {
    ['dtrap-f','dtrap-a','dtrap-b','dtrap-c','dtrap-d'].forEach(id=>{const el=document.getElementById(id);if(el){el.value='';el.classList.remove('has-error');}});
    ['dtrap-m','dtrap-n'].forEach(id=>{const el=document.getElementById(id);if(el)el.value='4';});
    U.hide('dtrap-results'); U.clearErr('dtrap-error');
    const btn=document.getElementById('dt-toggle-btn');
    if(btn) btn.style.display='none';
    this._gridData=null;
  },

  example() {
    document.getElementById('dtrap-f').value='x + y';
    document.getElementById('dtrap-a').value='0';
    document.getElementById('dtrap-b').value='1';
    document.getElementById('dtrap-c').value='0';
    document.getElementById('dtrap-d').value='1';
    document.getElementById('dtrap-m').value='4';
    document.getElementById('dtrap-n').value='4';
  }
};