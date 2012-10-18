%w[python-dev make gcc build-essential git-core python-pip].each do |p|
 package p
end

execute "pip" do
    command "sudo pip install nosetests"
    action :run
end
