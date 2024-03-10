template<class Type>
struct vmint {
  vector<matrix<int> > x;
  vmint(SEXP r_list) {
    x.resize(LENGTH(r_list));
    for(int i = 0; i < x.size(); i++) {
      x(i) = asMatrix<int>(VECTOR_ELT(r_list, i));
    }
  }
  int size() { return x.size(); }
  matrix<int> operator() (int i) { return x(i); }
};

template<class Type>
struct vvint {
  vector<vector<int> > x;
  vvint(SEXP r_list) {
    x.resize(LENGTH(r_list));
    for(int i = 0; i < x.size(); i++) {
      x(i) = asVector<int>(VECTOR_ELT(r_list, i));
    }
  }
  int size() { return x.size(); }
  vector<int> operator() (int i) { return x(i); }
};
