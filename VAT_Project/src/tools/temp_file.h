
#ifndef TEMP_FILE_H_
#define TEMP_FILE_H_

class TempFile
{
	public:
	TempFile():
		file_des_ (-1)
	{
		char *name = new char[VATParameters::tmpdir.length()+20];
		sprintf(name, "%s/VAT-tmp-XXXXXX", VATParameters::tmpdir.c_str());
		file_des_ = mkstemp(name);
		if(file_des_ == -1)
			throw std::runtime_error("Error opening temporary file.");
		unlink(name);
		delete[] name;
	}

	int file_descriptor() const
	{ return file_des_; }

private:

	int file_des_;

};

#endif /* TEMP_FILE_H_ */
