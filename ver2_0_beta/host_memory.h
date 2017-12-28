#pragma once
#ifndef HOST_MEMORY_H_
#define HOST_MEMORY_H_

/*------------------------------------------------------
* INCLUDES
*------------------------------------------------------*/
#include <stdlib.h>
#include <memory>


/*------------------------------------------------------
* FUNCTIONS
*------------------------------------------------------*/
template <typename data_type>
data_type* host_malloc(size_t mem_size)
{
	return (data_type*)malloc(mem_size * sizeof(data_type));
}

template<typename Derived, typename Base, typename Del>
std::unique_ptr<Derived, Del> static_unique_ptr_cast(std::unique_ptr<Base, Del>&& p)
{
	auto d = static_cast<Derived *>(p.release());
	return std::unique_ptr<Derived, Del>(d, std::move(p.get_deleter()));
}

template<typename Derived, typename Base, typename Del>
std::unique_ptr<Derived, Del> dynamic_unique_ptr_cast(std::unique_ptr<Base, Del>&& p)
{
	if (Derived *result = dynamic_cast<Derived *>(p.get())) {
		p.release();
		return std::unique_ptr<Derived, Del>(result, std::move(p.get_deleter()));
	}
	return std::unique_ptr<Derived, Del>(nullptr, p.get_deleter());
}

template<typename Derived, typename Base>
std::unique_ptr<Derived> dynamic_unique_ptr_cast(std::unique_ptr<Base>&& p)
{
	if (Derived *result = dynamic_cast<Derived *>(p.get())) {
		p.release();
		return std::unique_ptr<Derived>(result);
	}
	return std::unique_ptr<Derived>(nullptr);
}

#endif /* !HOST_MEMORY_H_*/
