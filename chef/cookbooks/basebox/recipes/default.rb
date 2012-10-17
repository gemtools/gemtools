%w[vim python-dev make gcc build-essential git-core python-pip samtools].each do |p|
 package p
end

execute "pip" do
    command "sudo pip install nose"
    action :run
end
