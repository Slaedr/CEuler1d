#ifndef __ARRAY_H
#define __ARRAY_H

#include <cstdio>

template<typename T>
class Array1d
{
protected:
	T* data;
	size_t size;
	bool alloc;
	bool device_alloc;
public:
	Array1d()
	{
		size = 0;
		alloc = false;
		device_alloc = false;
	}
	Array1d(size_t sz)
	{
		size = sz;
		data = new T[size];
		alloc = true;
		device_alloc = false;
	}
	~Array1d()
	{
		if(device_alloc)
		{
			#pragma acc exit data delete(data[:size])
			#pragma acc exit data delete(this)
		}
		if(alloc)
			delete [] data;
	}
	
	void allocate(const size_t sz)
	{
		size = sz;
		data = new T[size];
		alloc = true;
	}

	void device_allocate()
	{
		if(alloc)
		{
			#pragma acc enter data copyin(this)
			#pragma acc enter data create(data[:size])
			device_alloc = true;
		}
	}

	void device_delete()
	{
		if(device_alloc)
		{
			#pragma acc exit data delete(data[:size])
			#pragma acc exit data delete(this)
			device_alloc = false;
		}
	}

	void device_update()
	{
		if(alloc && device_alloc)
		{
			#pragma acc update device(data[:size])
		}
	}

	void update()
	{
		if(alloc && device_alloc)
		{
			#pragma acc update self(data[:size])
		}
	}

#pragma acc routine seq
	T& operator[](const size_t i)
	{
		return data[i];
	}

#pragma acc routine seq
	T operator[](const size_t i) const
	{
		return data[i];
	}
};

/** \brief Row-major 2D array storage
 *
 * Stores members in a contiguous 1D block.
 */
template <typename T>
class Array2d
{
protected:
	T* data;
	size_t nrows;
	size_t ncols;
	size_t size;
	bool alloc;
public:
	Array2d()
	{
		nrows = ncols = size = 0;
		alloc = false;
	}
	Array2d(const size_t nr, const size_t nc)
	{
		nrows = nr; ncols = nc;
		size = ncols*nrows;
		data = new T[size];
		alloc = true;
	}
	void allocate(const size_t nr, const size_t nc)
	{
		nrows = nr; ncols = nc;
		size = ncols*nrows;
		data = new T[size];
		alloc = true;
	}
	~Array2d()
	{
		if(alloc)
			delete [] data;
	}

	T& operator[][](const size_t i, const size_t j)
	{
		return data[i*nrows+j];
	}

	T operator[][](const size_t i, const size_t j) const
	{
		return data[i*nrows+j];
	}
};

#endif
